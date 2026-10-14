/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Held loci follow their set (seqspec:interaction#selection): on every assembled set, the selection, the highlight
// channels and an external hover are carried to it. Edge-rate: once per held loci per new set.

import { Loci, type SeqSet } from "@molstar/seq-core/model";

/** Any loci the viewer holds. */
export type SeqLociAny = Loci;

/** Within a lineage, whether every owner entity of `l` is the same in `b` as in `a` (its indices then read the same). */
function sameOwners(l: Loci, a: SeqSet, b: SeqSet): boolean {
  if ("alignment" in l) return a.alignments[l.alignment] === b.alignments[l.alignment];
  for (const e of l.elements) {
    if ("annotation" in e) {
      if (a.annotations[e.annotation] !== b.annotations[e.annotation]) return false;
      continue;
    }
    const u = e.unit;
    const [ua, ub] = [a.units, b.units];
    if (!ub.alive[u] || a.molecules[ua.molecule[u]] !== b.molecules[ub.molecule[u]]) return false;
    if (ua.rangeStart[u] !== ub.rangeStart[u] || ua.rangeStop[u] !== ub.rangeStop[u]) return false;
  }
  return true;
}

/**
 * `loci` in `set`: itself when it is already there; within a lineage, rebound to `set` when its owners are the same
 * objects, else remapped by key (`Loci.remap`: an edited range or item table read anew, tombstoned ids dropped);
 * across lineages remapped by key. Null: nothing of it is left.
 */
export function carry<L extends Loci>(loci: L, set: SeqSet): L | null {
  const from = loci.set;
  if (from === set) return loci;
  if (from.lineage === set.lineage && sameOwners(loci, from, set)) return { ...loci, set };
  const r = Loci.remap(loci, from, set);
  return Loci.isEmpty(r) ? null : r;
}
