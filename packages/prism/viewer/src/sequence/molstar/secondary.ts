/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `molstar:secondary_structure` (seqspec:libraries#the-mol-adapter): a structure's helices, strands and turns as a
// `ranges` annotation on its complex's polymer units. The assignment the host's 3D view draws (Mol*'s computed secondary
// structure on the drawn Structure: the file's where it has one, else DSSP, as the cartoon shows it), else the one the
// model carries (mmCIF's `struct_conf` and `struct_sheet_range`). Coil is no range. Nothing the 3D view would not
// compute itself, and no network.

import {
  nullable,
  OptionalField,
  type ProviderEntry,
  RequiredField,
  ref,
  SimpleParamsSchema,
  str,
} from '@molstar/seq-core/format';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { ModelSecondaryStructure } from '@molstar/model/formats/structure/property/secondary-structure';
import { SecondaryStructureProvider as ComputedSecondaryStructure } from '@molstar/model/props/computed/secondary-structure';
import type { AtomicHierarchy } from '@molstar/model/model/structure/model/properties/atomic';
import type { SecondaryStructure } from '@molstar/model/model/structure/model/properties/secondary-structure';
import { type Structure, Unit } from '@molstar/model/model/structure';
import type { ResidueIndex } from '@molstar/model/model/structure/model/indexing';
import { SecondaryStructureType } from '@molstar/model/model/structure/model/types';
import { chainFor, complexUnits } from './seq-units.js';
import { StructureSource, structureOf } from './source.js';

/** A range's `category` in `molstar:secondary_structure`'s annotation. */
export type SecondaryCategory = 'helix' | 'strand' | 'turn';

/**
 * A residue's class from Mol*'s flags: a helix (any), a strand (beta), a turn or bend; null: coil. `bends` false: a bend
 * is coil (a computed assignment's: DSSP assigns them freely, and the cartoon draws none; a file lists the ones it means).
 */
export function secondaryCategory(type: SecondaryStructureType, bends = true): SecondaryCategory | null {
  if (
    SecondaryStructureType.is(type, SecondaryStructureType.Flag.Helix) ||
    SecondaryStructureType.is(type, SecondaryStructureType.Flag.DoubleHelix)
  )
    return 'helix';
  if (SecondaryStructureType.is(type, SecondaryStructureType.Flag.Beta)) return 'strand';
  if (
    SecondaryStructureType.is(type, SecondaryStructureType.Flag.Turn) ||
    (bends && SecondaryStructureType.is(type, SecondaryStructureType.Flag.Bend))
  )
    return 'turn';
  return null;
}

/**
 * One row per element of the assignment (a helix, strand, turn or a file's bend) on each polymer unit of the complex `of`: two elements
 * that touch stay two rows (a kinked helix given as two helices is two). The residues are the source's first model as
 * deposited (the asymmetric unit), as `molstar:plddt` reads them, whatever assembly the complex shows; the assignment is
 * per `label_asym_id`, so an assembly's copy of a chain gets its chain's rows. It is the drawn Structure's computed one
 * where the drawn chain has the same residues (its model's hierarchy: another model of a trajectory may not), else the
 * model's own (`struct_conf`, `struct_sheet_range`). A chain with none has no rows: the annotation can be empty, never an
 * error.
 *
 * `copies: "one"` (a section that merges its copies): one row set per molecule, on the molecule (it reaches every copy),
 * merged from its chains' so a merged block shows one, not every chain's stacked (3J3Q: 1,356 chains, each its own DSSP).
 * Per residue the class most of its chains observing it have, the first chain's on a tie; an element ends where most of
 * the chains with its class start a new one.
 */
export const SecondaryStructureProvider: ProviderEntry = {
  name: 'molstar:secondary_structure',
  options: SimpleParamsSchema({
    of: RequiredField(ref('complex'), 'the complex of a structure'),
    // Not read: the host writes the drawn assignment's settings here, so a change of them recomputes the annotation.
    assignment: OptionalField(nullable(str), null, "the 3D view's assignment settings, as the host knows them"),
    copies: OptionalField(
      nullable(str),
      null,
      '"one": each molecule\'s chains merged into one row set on the molecule (a section merging copies); else per chain',
    ),
  }),
  geometry: 'ranges',
  async compute(ctx, options) {
    const { of, copies } = options as { of: string; copies?: string | null };
    const parsed = ctx.source(of);
    const src = parsed?.kind === 'object' ? parsed.value : null;
    const set = ctx.set;
    if (!StructureSource.is(src) || !set) return null;
    const structure = await structureOf(src, 0, null);
    const drawn = await drawnAssignment(src.drawn);
    const target: string[] = [];
    const start: number[] = [];
    const stop: number[] = [];
    const category: string[] = [];
    const merged = copies === 'one' ? new Map<number, Votes>() : null;
    for (const u of complexUnits(set, of)) {
      const c = chainFor(structure, u.mvs);
      if (c?.kind !== 'polymer') continue;
      // The drawn assignment indexes the drawn model's residues: only where they are this model's (the same hierarchy).
      const d = drawn?.get(c.labelAsymId);
      const own = ModelSecondaryStructure.Provider.get(c.unit.model);
      const ss = d && d.hierarchy === c.unit.model.atomicHierarchy ? d.ss : own;
      if (!ss) continue;
      if (!merged) {
        elements(c, ss, ss === own, u.first, null, (kind, from, to) => {
          target.push(u.ref);
          start.push(from);
          stop.push(to);
          category.push(Categories[kind]);
        });
        continue;
      }
      const mol = set.units.molecule[u.unit];
      let v = merged.get(mol);
      if (!v) {
        v = votesFor(set.molecules[mol]?.length ?? 0);
        merged.set(mol, v);
      }
      const votes = v;
      const firstChain = votes.chains++ === 0;
      elements(
        c,
        ss,
        ss === own,
        u.first,
        (at) => {
          if (at < 0 || at >= votes.length) return;
          votes.observed[at]++;
          if (firstChain) votes.first[at] = 0;
        },
        (kind, from, to) => {
          for (let at = Math.max(0, from); at < Math.min(to, votes.length); at++) {
            votes.votes[3 * at + kind]++;
            if (firstChain) votes.first[at] = kind + 1;
          }
          if (from >= 0 && from < votes.length) votes.starts[3 * from + kind]++;
        },
      );
    }
    for (const [mol, v] of merged ?? []) {
      const key = set.strings.get(set.molecules[mol]?.key ?? -1);
      consensus(v, (kind, from, to) => {
        target.push(key);
        start.push(from);
        stop.push(to);
        category.push(Categories[kind]);
      });
    }
    return { key: '', geometry: 'ranges', label: 'Secondary structure', columns: { target, start, stop, category } };
  },
};

/** The classes rows are of, by index (`elements`, `Votes`). */
const Categories: readonly SecondaryCategory[] = ['helix', 'strand', 'turn'];
const categoryIndex = (c: SecondaryCategory) => (c === 'helix' ? 0 : c === 'strand' ? 1 : 2);

/**
 * A chain's elements as `element(class, from, to)`, positions on its molecule (label_seq_id less the numbering's
 * first): a run goes on while its element does and the residues are consecutive; `residue(at)` per residue observed.
 */
function elements(
  c: { count: number; residues: ArrayLike<number>; labels: ArrayLike<number> },
  ss: SecondaryStructure,
  fileOwn: boolean,
  first: number,
  residue: ((at: number) => void) | null,
  element: (kind: number, from: number, to: number) => void,
) {
  let open: SecondaryCategory | null = null;
  let key = -1;
  let from = 0;
  let last = 0;
  const close = () => {
    if (open === null) return;
    element(categoryIndex(open), from, last + 1);
    open = null;
  };
  for (let x = 0; x < c.count; x++) {
    // The assignment's own index of the residue (`getIndex`), and its element there (`key`).
    const i = ss.getIndex(c.residues[x] as ResidueIndex);
    // A file's bends are rows; DSSP's are not (Mol*'s computed assignment is the file's own object where it has one).
    const kind = secondaryCategory(ss.type[i], fileOwn);
    const k = ss.key[i];
    const at = c.labels[x] - first;
    residue?.(at);
    if (kind !== open || k !== key || (open !== null && at !== last + 1)) {
      close();
      if (kind !== null) {
        open = kind;
        key = k;
        from = at;
      }
    }
    last = at;
  }
  close();
}

/**
 * A molecule's chains' assignments counted per residue: chains observing it, chains with each class there and chains
 * whose element of that class starts there (3 per residue, `Categories`' order), and the first chain's (0: coil, else
 * its class + 1; -1: not observed).
 */
interface Votes {
  readonly length: number;
  chains: number;
  readonly observed: Int32Array;
  readonly votes: Int32Array;
  readonly starts: Int32Array;
  readonly first: Int8Array;
}
const votesFor = (length: number): Votes => ({
  length,
  chains: 0,
  observed: new Int32Array(length),
  votes: new Int32Array(3 * length),
  starts: new Int32Array(3 * length),
  first: new Int8Array(length).fill(-1),
});

/** The chains' merged elements (`SecondaryStructureProvider`, `copies: "one"`). */
function consensus(v: Votes, element: (kind: number, from: number, to: number) => void) {
  let open = -1;
  let from = 0;
  for (let at = 0; at <= v.length; at++) {
    // The class most chains observing the residue have (coil: -1), the first chain's on a tie, else the lowest.
    let best = -1;
    if (at < v.length && v.observed[at] > 0) {
      const count = (k: number) =>
        k < 0 ? v.observed[at] - v.votes[3 * at] - v.votes[3 * at + 1] - v.votes[3 * at + 2] : v.votes[3 * at + k];
      let most = count(-1);
      for (let k = 0; k < 3; k++) most = Math.max(most, count(k));
      const mine = v.first[at] - 1;
      if (mine >= -1 && count(mine) === most) best = mine;
      else while (count(best) !== most) best++;
    }
    const split = best >= 0 && best === open && 2 * v.starts[3 * at + best] > v.votes[3 * at + best];
    if (best !== open || split) {
      if (open >= 0) element(open, from, at);
      open = best;
      from = at;
    }
  }
}

/**
 * The secondary structure the 3D view draws, by `label_asym_id`, with the hierarchy its residue indices are of: Mol*'s
 * computed property on the drawn Structure (per unit symmetry group), attached with its props as they are when the 3D
 * view has not yet (the same computation its cartoon makes, done once: attaching again with equal props does nothing).
 * Undefined: no drawn Structure.
 */
async function drawnAssignment(
  drawn: Structure | undefined,
): Promise<Map<string, { ss: SecondaryStructure; hierarchy: AtomicHierarchy }> | undefined> {
  if (!drawn) return undefined;
  const root = drawn.root;
  if (!ComputedSecondaryStructure.get(root).value) {
    await ComputedSecondaryStructure.attach(
      { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() },
      root,
    );
  }
  const value = ComputedSecondaryStructure.get(root).value;
  if (!value) return undefined;
  const byChain = new Map<string, { ss: SecondaryStructure; hierarchy: AtomicHierarchy }>();
  for (const group of root.unitSymmetryGroups) {
    const u = group.units[0];
    if (!Unit.isAtomic(u)) continue;
    const ss = value.get(u.invariantId);
    if (!ss) continue;
    const h = u.model.atomicHierarchy;
    const asym = h.chains.label_asym_id.value(h.chainAtomSegments.index[u.elements[0]]);
    if (!byChain.has(asym)) byChain.set(asym, { ss, hierarchy: h });
  }
  return byChain;
}
