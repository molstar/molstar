/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Translator's Mol* side (seqspec:libraries#the-mol-adapter; seqspec:libraries#translate item 4): a loaded
// structure registered as an external complex system, `molstar:<ref>`, whose members are its chain instances,
// `molstar:<ref>/<label_asym_id>@<instance_id>`, each indexed by element (a polymer's residues, any other chain's heavy
// atoms) and mapped onto the Seq* unit of the same chain: by label_seq_id for a polymer, by heavy-atom order for the
// others. Mol* loci become runs in those systems and back; `translate` carries runs on to Seq*'s systems, or the viewer
// projects them as any external system's.

import {
  type CoordinateMap,
  type CoordinateSystem,
  CoordinateSystemKind,
  Extents,
  KeyKind,
  type SeqSet,
} from '@molstar/seq-core/model';
import { TranslateMaps, type Translator } from '@molstar/seq-core/translate';
import { OrderedSet } from 'molstar/lib/mol-data/int';
import { type Structure, StructureElement, type Unit } from 'molstar/lib/mol-model/structure';
import { type ChainInstance, chainsOf } from './chains.js';
import { complexUnits, seqRoutes } from './seq-units.js';

const unitsOf = (chains: readonly LinkedChain[]) => {
  const units = new Set<number>();
  for (const c of chains) if (c.unit) units.add(c.unit.id);
  return units;
};

export interface LinkedChain {
  readonly key: string;
  readonly system: CoordinateSystem;
  readonly chain: ChainInstance;
  /** The Seq* unit it maps onto (`{ kind: unit, id }`); null: none in the complex. */
  readonly unit: CoordinateSystem | null;
  readonly unitRef: string | null;
}

const UNIT = CoordinateSystemKind.unit;

const sameRun = (e: Extents, a: number, b: number) =>
  e.kind[a] === e.kind[b] &&
  e.id[a] === e.id[b] &&
  e.row[a] === e.row[b] &&
  e.start[a] === e.start[b] &&
  e.stop[a] === e.stop[b];

/** One structure linked to one Seq* complex through a `Translator` (the viewer's `ctx.systems`, or a host's own). */
export class StructureLink {
  readonly key: string;
  readonly system: CoordinateSystem;
  readonly chains: readonly LinkedChain[];
  readonly structure: Structure;
  /** The set its units are in: a query's runs name them (query.ts). */
  readonly set: SeqSet;
  private readonly t: Translator;
  /** Member index by external system id. */
  private readonly byId = new Map<number, number>();
  /** Mol* unit id → chain instance per position, element per position (chainsOf's). */
  private readonly chainAt: ReadonlyMap<number, Int32Array>;
  private readonly elementAt: ReadonlyMap<number, Int32Array>;
  /** Chain instance (chainsOf's index) → member, -1: not linked (water). */
  private readonly chainMember: Int32Array;
  private readonly one = Extents.create(1);
  private readonly scratch = Extents.create(16);

  private constructor(t: Translator, ref: string, structure: Structure, set: SeqSet, complex: string) {
    this.t = t;
    this.structure = structure;
    this.set = set;
    this.key = `molstar:${ref}`;
    const index = chainsOf(structure);
    this.chainAt = index.chainAt;
    this.elementAt = index.elementAt;
    this.chainMember = new Int32Array(index.chains.length).fill(-1);
    const seqUnits = complexUnits(set, complex);
    const target = { kind: 0, id: 0, row: -1 };
    const chains: LinkedChain[] = [];
    index.chains.forEach((c, k) => {
      if (c.kind === 'water') return;
      const u =
        seqUnits.find((s) => s.mvs.label_asym_id === c.labelAsymId && s.mvs.instance_id === c.instanceId) ??
        seqUnits.find((s) => s.mvs.label_asym_id === c.labelAsymId);
      const unit =
        u && set.byName(u.ref, target) && target.kind === KeyKind.unit ? { kind: UNIT, id: target.id, row: -1 } : null;
      const key = `${this.key}/${c.labelAsymId}@${c.instanceId}`;
      const maps = (s: CoordinateSystem): CoordinateMap[] => {
        if (!unit || !u) return [];
        const mol = set.molecules[set.units.molecule[unit.id]];
        const length = mol?.length ?? 0;
        if (c.kind !== 'polymer')
          return [TranslateMaps.bySegments(s, unit, Int32Array.of(0, 0, Math.min(c.count, length)))];
        const labels = new Int32Array(length);
        for (let i = 0; i < length; i++) labels[i] = i + u.first;
        return [TranslateMaps.byLabels(s, c.labels, unit, labels)];
      };
      const system = t.register(
        {
          key,
          kind: 'external',
          start: 0,
          stop: c.count,
          refs: u ? { unit: u.ref, complex } : {},
          members: null,
          structure: {
            structureRef: ref,
            entity: c.entityId,
            labelAsymId: c.labelAsymId,
            authAsymId: c.authAsymId,
            instanceId: c.instanceId,
          },
          numberings: [],
          owner: 'molstar',
        },
        maps,
      );
      this.chainMember[k] = chains.length;
      this.byId.set(system.id, chains.length);
      chains.push({ key, system, chain: c, unit, unitRef: u?.ref ?? null });
    });
    this.chains = chains;
    this.system = t.register(
      {
        key: this.key,
        kind: 'external',
        start: 0,
        stop: 0,
        refs: { complex },
        members: chains.map((c) => c.key),
        structure: { structureRef: ref },
        numberings: [],
        owner: 'molstar',
      },
      seqRoutes(set, unitsOf(chains)),
    );
  }

  /** Registers `structure`'s chains as `molstar:<ref>` and its members, mapped onto `complex`'s units in `set`. */
  static attach(t: Translator, ref: string, structure: Structure, set: SeqSet, complex: string): StructureLink {
    return new StructureLink(t, ref, structure, set, complex);
  }

  /** Unregisters the systems; their keys keep their ids, so a later `attach` under `ref` names the same systems. */
  detach(): void {
    for (const c of this.chains) this.t.unregister(c.key);
    this.t.unregister(this.key);
  }

  /**
   * A Mol* loci (this structure's) as runs in the member systems, appended to `out`: consecutive elements of one member
   * in one run. Allocation-free once `out` is large enough (the 3D hover's rate).
   */
  fromLoci(loci: StructureElement.Loci, out: Extents): Extents {
    for (let k = 0; k < loci.elements.length; k++) {
      const { unit, indices } = loci.elements[k];
      const cAt = this.chainAt.get(unit.id);
      const eAt = this.elementAt.get(unit.id);
      if (!cAt || !eAt) continue;
      for (let i = 0, n = OrderedSet.size(indices); i < n; i++) {
        const p = OrderedSet.getAt(indices, i);
        const c = cAt[p];
        const m = c < 0 ? -1 : this.chainMember[c];
        if (m < 0) continue;
        const e = eAt[p];
        const s = this.chains[m].system;
        const last = out.count - 1;
        const same = last >= 0 && out.kind[last] === s.kind && out.id[last] === s.id;
        if (same && e >= out.start[last] && e < out.stop[last]) continue;
        if (same && e === out.stop[last]) out.stop[last] = e + 1;
        else Extents.add(out, s, e, e + 1);
      }
    }
    return out;
  }

  /**
   * Runs in any system with a route to this structure (a Seq* unit's, a molecule's, a member's) as runs in the member
   * systems, appended to `out`.
   */
  toMembers(runs: Extents, out: Extents): Extents {
    for (let i = 0; i < runs.count; i++) {
      if (runs.kind[i] === CoordinateSystemKind.external && this.byId.has(runs.id[i])) {
        Extents.add(out, this.chains[this.byId.get(runs.id[i]) as number].system, runs.start[i], runs.stop[i]);
        continue;
      }
      const from = { kind: runs.kind[i], id: runs.id[i], row: runs.row[i] };
      this.one.count = 0;
      Extents.add(this.one, from, runs.start[i], runs.stop[i]);
      this.t.translate(this.one, from, this.system, out);
    }
    return out;
  }

  /** The member a run of `toMembers`' names; null: not this structure's. */
  memberOf(runs: Extents, i: number): LinkedChain | null {
    const m = runs.kind[i] === CoordinateSystemKind.external ? this.byId.get(runs.id[i]) : undefined;
    return m === undefined ? null : this.chains[m];
  }

  /** Runs in any system with a route to this structure as one Mol* loci; edge rate, allocates. */
  toLoci(runs: Extents): StructureElement.Loci {
    const members = this.scratch;
    members.count = 0;
    this.toMembers(runs, members);
    const byUnit = new Map<number, { unit: Unit; at: number[] }>();
    for (let i = 0; i < members.count; i++) {
      const c = this.memberOf(members, i)?.chain;
      if (!c) continue;
      const entry = byUnit.get(c.unit.id) ?? { unit: c.unit, at: [] };
      byUnit.set(c.unit.id, entry);
      for (let e = Math.max(0, members.start[i]); e < Math.min(c.count, members.stop[i]); e++) {
        for (let p = c.start[e]; p < c.end[e]; p++) entry.at.push(p);
      }
    }
    const elements: { unit: Unit; indices: OrderedSet<StructureElement.UnitIndex> }[] = [];
    for (const unit of this.structure.units) {
      const entry = byUnit.get(unit.id);
      if (!entry || entry.at.length === 0) continue;
      const at = Int32Array.from(new Set(entry.at)).sort();
      const contiguous = at[at.length - 1] - at[0] + 1 === at.length;
      const indices = contiguous ? OrderedSet.ofBounds(at[0], at[at.length - 1] + 1) : OrderedSet.ofSortedArray(at);
      elements.push({ unit, indices: indices as OrderedSet<StructureElement.UnitIndex> });
    }
    return StructureElement.Loci(this.structure, elements);
  }

  /** Member runs as runs in their Seq* units, appended to `out` (the `structure` input of a SeqQL query). */
  toSeq(runs: Extents, out: Extents): Extents {
    for (let i = 0; i < runs.count; i++) {
      const linked = this.memberOf(runs, i);
      if (!linked?.unit) continue;
      this.one.count = 0;
      Extents.add(this.one, linked.system, runs.start[i], runs.stop[i]);
      const before = out.count;
      this.t.translate(this.one, linked.system, linked.unit, out);
      // Copies of a chain (an assembly's) land on one unit: a run equal to the one before it is dropped.
      if (out.count === before + 1 && before > 0 && sameRun(out, before, before - 1)) out.count--;
    }
    return out;
  }

  /** A Mol* loci as runs in Seq* units: `fromLoci`, then `toSeq`. */
  seqRunsOf(loci: StructureElement.Loci, out: Extents): Extents {
    const members = this.scratch;
    members.count = 0;
    return this.toSeq(this.fromLoci(loci, members), out);
  }
}
