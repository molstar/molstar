/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Coordinate systems at the edge (seqspec:interaction#focus; seqspec:libraries#translate): `ctx.systems`, one
// `Translator` holding a host's systems and the maps registered with them, which also names the shown set's own systems
// by key, `seqstar:<ref>` (a molecule, a unit, a complex, an alignment's columns, a row `<alignment>/<k>`). Their
// descriptors are built on first use, once per system and lineage, and cached by key, so a request naming one reads a
// cache; `coordinatesOf` gives loci as `Coordinates` for events. Edge rate but the reads.

import {
  type CoordinateMap,
  type CoordinateSystem,
  Extents,
  CoordinateSystemKind as K,
  KeyKind,
  type Loci,
  type SeqSet,
  Systems,
} from '@molstar/seq-core/model';
import { type SystemDescriptor, Translator } from '@molstar/seq-core/translate';
import type { Coordinates } from './projection.js';
import { state } from './state.js';

const PREFIX = 'seqstar:';
/** An own system's descriptor, valid while `owner` (the table or entity it reads) is the set's. */
interface Own {
  readonly set: SeqSet;
  readonly owner: unknown;
  readonly system: CoordinateSystem;
  readonly descriptor: SystemDescriptor;
}
/** The ref naming a Seq* system in `set` (core's key, `msa` for its columns); null: a kind no key names. */
const refOf = Systems.refOf;
const domain = new Int32Array(2);

export class SystemRegistry extends Translator {
  /** Bumped by every registration: channels holding `Coordinates` light them again. */
  readonly version = state(0);
  /** Own descriptors built so far (a counter for tests). */
  built = 0;
  private set: SeqSet | null = null;
  private readonly own = new Map<string, Own>();
  private readonly key = { kind: 0, id: 0, row: -1 };

  /** The shown set (`graph.assembled`): own descriptors are kept within its lineage while what they read is kept. */
  update(set: SeqSet | null): void {
    if (!set || set === this.set) return;
    if (this.set?.lineage !== set.lineage) this.own.clear();
    this.set = set;
  }

  override register(
    d: SystemDescriptor,
    maps: readonly CoordinateMap[] | ((s: CoordinateSystem) => readonly CoordinateMap[]) = [],
    system?: CoordinateSystem,
  ): CoordinateSystem {
    const s = super.register(d, maps, system);
    this.version.set(this.version.value + 1);
    return s;
  }
  override unregister(key: string): void {
    super.unregister(key);
    this.version.set(this.version.value + 1);
  }
  override addMap(map: CoordinateMap): void {
    super.addMap(map);
    this.version.set(this.version.value + 1);
  }

  /** A host's key, else the shown set's own (`seqstar:<ref>`). */
  override system(key: string): CoordinateSystem | null {
    return super.system(key) ?? this.ownOf(key)?.system ?? null;
  }
  override descriptor(key: string): SystemDescriptor | null {
    return super.descriptor(key) ?? this.ownOf(key)?.descriptor ?? null;
  }
  override descriptorOf(system: CoordinateSystem): SystemDescriptor | null {
    const d = super.descriptorOf(system);
    const ref = d || !this.set ? null : refOf(this.set, system.kind, system.id, system.row);
    return d ?? (ref === null ? null : (this.ownOf(PREFIX + ref)?.descriptor ?? null));
  }

  /**
   * Loci as `Coordinates` in the shown set's systems (edge rate, allocates): each element's runs under its system's
   * key, an item as its rows' ranges; a unit of a complex also lists the complex, whose descriptor's `members` name its
   * units, so a host builds one structure loci from the runs.
   */
  coordinatesOf(loci: Loci): Coordinates {
    const set = loci.set;
    const ext = Extents.create(8);
    Extents.addLoci(ext, loci);
    const systems: string[] = [];
    const runs: number[] = [];
    const index = (ref: string) => {
      const k = systems.indexOf(PREFIX + ref);
      return k >= 0 ? k : systems.push(PREFIX + ref) - 1;
    };
    for (let i = 0; i < ext.count; i++) {
      const ref = refOf(set, ext.kind[i], ext.id[i], ext.row[i]);
      if (ref !== null) runs.push(index(ref), ext.start[i], ext.stop[i]);
      const cx = ext.kind[i] === K.unit ? set.complexes[set.units.complex[ext.id[i]]] : null;
      if (cx) index(set.strings.get(cx.key));
    }
    return { systems, runs: Int32Array.from(runs) };
  }

  /** An own key's system and descriptor, cached; null: not `seqstar:`, or naming nothing in the shown set. */
  private ownOf(key: string): Own | null {
    const set = this.set;
    const hit = this.own.get(key);
    if (!set || (hit && hit.set === set)) return hit ?? null;
    const t = this.key;
    if (!key.startsWith(PREFIX) || !set.byName(key.slice(PREFIX.length), t)) return null;
    // By `KeyKind` (molecule, unit, alignment, row, mapping, annotation, complex): the coordinate system it names.
    const kind = [K.molecule, K.unit, K.columns, K.row, -1, -1, K.complex][t.kind];
    const id = t.id;
    const owner =
      t.kind === KeyKind.unit
        ? set.units
        : t.kind === KeyKind.molecule
          ? set.molecules[id]
          : t.kind === KeyKind.complex
            ? set.complexes[id]
            : set.alignments[id];
    if (kind < 0 || !owner) return null;
    if (hit?.owner === owner) {
      this.own.set(key, { ...hit, set });
      return this.own.get(key) ?? null;
    }
    const own = { set, owner, system: { kind, id, row: t.row }, descriptor: this.describe(set, key, kind, id, t.row) };
    this.built++;
    this.own.set(key, own);
    return own;
  }

  private describe(set: SeqSet, key: string, kind: number, id: number, row: number): SystemDescriptor {
    const u = set.units;
    const str = (k: number) => set.strings.get(k);
    Systems.domainOf(set, kind, id, row, domain);
    const [start, stop] = [domain[0], Math.max(domain[0], domain[1])];
    const base = { key, start, stop, refs: {}, members: null, structure: null, numberings: [], owner: 'seqstar' };
    if (kind === K.unit) {
      const m = set.molecules[u.molecule[id]];
      const cx = set.complexes[u.complex[id]];
      const mvs = (u.custom[id]?.doc as { mvs?: SystemDescriptor['structure'] } | undefined)?.mvs ?? null;
      return {
        ...base,
        kind: 'unit',
        refs: { unit: str(u.key[id]), molecule: m ? str(m.key) : undefined, complex: cx ? str(cx.key) : undefined },
        structure: mvs,
        numberings: Systems.numberingNames(set, kind, id, row),
      };
    }
    if (kind === K.molecule) return { ...base, kind: 'molecule', refs: { molecule: key.slice(PREFIX.length) } };
    if (kind === K.complex) {
      const cx = set.complexes[id];
      const members = Array.from(cx?.units ?? [], (x) => PREFIX + str(u.key[x]));
      return { ...base, kind: 'complex', refs: { complex: key.slice(PREFIX.length) }, members };
    }
    const a = set.alignments[id];
    const alignment = a ? str(a.key) : undefined;
    if (kind === K.columns) return { ...base, kind: 'columns', refs: { alignment } };
    return { ...base, kind: 'row', refs: { alignment, row } };
  }
}
