/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A structure's chain instances (seqspec:libraries#the-mol-adapter): one per label_asym_id and operator, whatever
// Mol*'s units are (Mol* merges waters and consecutive one-atom chains into one unit). Each lists its elements as Seq*
// indexes them: residues for a polymer chain, heavy atoms for any other. Built once per structure, read by the complex
// provider, the structure providers and the link.

import { type ElementIndex, type Structure, Unit } from '@molstar/model/model/structure';

export type ChainKind = 'polymer' | 'ligand' | 'branched' | 'water';

export interface ChainInstance {
  readonly labelAsymId: string;
  readonly authAsymId: string;
  readonly entityId: string;
  /** Mol*'s operator `instanceId` (MVS's `instance_id`): "1_555" in a deposited model, "ASM-1"-like in an assembly. */
  readonly instanceId: string;
  readonly identity: boolean;
  readonly kind: ChainKind;
  /** The Mol* unit holding it. */
  readonly unit: Unit.Atomic;
  readonly count: number;
  /** Per element, its atoms' positions in `unit.elements`: [start, end). */
  readonly start: Int32Array;
  readonly end: Int32Array;
  /** Per element: a polymer's label_seq_id; else its ordinal + 1. */
  readonly labels: Int32Array;
  /** Per element, its residue. */
  readonly residues: Int32Array;
  /** All the chain's atoms in `unit.elements` (hydrogens and repeated names too): [start, end). */
  readonly span: readonly [start: number, end: number];
}

/** Per Mol* unit (by `unit.id`): each position's chain instance (an index into the list) and element, -1 for none. */
export interface ChainIndex {
  readonly chains: readonly ChainInstance[];
  readonly chainAt: ReadonlyMap<number, Int32Array>;
  readonly elementAt: ReadonlyMap<number, Int32Array>;
}

const cache = new WeakMap<Structure, ChainIndex>();
const KINDS: Record<string, ChainKind> = {
  polymer: 'polymer',
  'non-polymer': 'ligand',
  macrolide: 'ligand',
  branched: 'branched',
  water: 'water',
};

/** The structure's chain instances, in its units' order; atomic units only (coarse beads are not indexed). */
export function chainsOf(structure: Structure): ChainIndex {
  const hit = cache.get(structure);
  if (hit) return hit;
  const chains: ChainInstance[] = [];
  const chainAt = new Map<number, Int32Array>();
  const elementAt = new Map<number, Int32Array>();
  for (const unit of structure.units) {
    if (!Unit.isAtomic(unit)) continue;
    const els = unit.elements;
    const n = els.length;
    const cAt = new Int32Array(n).fill(-1);
    const eAt = new Int32Array(n).fill(-1);
    chainAt.set(unit.id, cAt);
    elementAt.set(unit.id, eAt);
    const h = unit.model.atomicHierarchy;
    const ent = unit.model.entities.data;
    const op = unit.conformation.operator;
    let p = 0;
    while (p < n) {
      const c = h.chainAtomSegments.index[els[p]];
      let q = p;
      while (q < n && h.chainAtomSegments.index[els[q]] === c) q++;
      const e = h.index.getEntityFromChain(c);
      const kind = KINDS[ent.type.value(e)] ?? 'ligand';
      const b = { start: [] as number[], end: [] as number[], labels: [] as number[], residues: [] as number[] };
      if (kind === 'polymer') {
        for (let i = p; i < q; ) {
          const r = h.residueAtomSegments.index[els[i]];
          let j = i;
          while (j < q && h.residueAtomSegments.index[els[j]] === r) j++;
          b.start.push(i);
          b.end.push(j);
          b.labels.push(h.residues.label_seq_id.value(r));
          b.residues.push(r);
          i = j;
        }
      } else if (kind === 'water') {
        // Listed with no elements: nothing links a water chain (the complex and the link skip them), and a Set per water
        // residue was a third of the index of a structure with 12k waters.
      } else {
        // Heavy atoms, one per name in a residue (an alternate location is not another atom).
        let seen = new Set<string>();
        let last = -1;
        for (let i = p; i < q; i++) {
          const a = els[i] as ElementIndex;
          const r = h.residueAtomSegments.index[a];
          if (r !== last) seen = new Set();
          last = r;
          const sym = h.atoms.type_symbol.value(a);
          const name = h.atoms.label_atom_id.value(a);
          if (sym === 'H' || sym === 'D' || seen.has(name)) continue;
          seen.add(name);
          b.start.push(i);
          b.end.push(i + 1);
          b.labels.push(b.labels.length + 1);
          b.residues.push(r);
        }
      }
      const k = chains.length;
      chains.push({
        labelAsymId: h.chains.label_asym_id.value(c),
        authAsymId: h.chains.auth_asym_id.value(c),
        entityId: ent.id.value(e),
        instanceId: op.instanceId,
        identity: op.isIdentity,
        kind,
        unit,
        count: b.start.length,
        start: Int32Array.from(b.start),
        end: Int32Array.from(b.end),
        labels: Int32Array.from(b.labels),
        residues: Int32Array.from(b.residues),
        span: [p, q],
      });
      for (let x = 0; x < b.start.length; x++) {
        for (let i = b.start[x]; i < b.end[x]; i++) {
          cAt[i] = k;
          eAt[i] = x;
        }
      }
      p = q;
    }
  }
  const index = { chains, chainAt, elementAt };
  cache.set(structure, index);
  return index;
}
