/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An augmented interval tree (seqspec:sources#levels-and-the-zoom-gate, "The index for scale"): intervals sorted by
// their low ends, and an implicit binary tree over that order whose nodes hold the largest high end below them. The
// intervals meeting
// [q[0], q[1]] are found in O((k + 1) log n) whatever their lengths, where a prefix maximum of high ends degenerates
// to a scan of every earlier interval once one long interval exists. Keys are doubles (pairs: block × 2^32 + index),
// read from typed arrays so none crosses a call. Prepare builds the tree; queries are allocation-free.

/** What `overlapping` and the end index visit (K4's counters): tree nodes entered, intervals reported, ends read. */
export const IndexVisits = { nodes: 0, intervals: 0, ends: 0 };

/** The tree over the high ends in sorted order: leaves from the first power of two ≥ n, each node its children's max. */
export function intervalTree(hiSorted: Float64Array): Float64Array {
  let size = 1;
  while (size < hiSorted.length) size *= 2;
  const t = new Float64Array(2 * size).fill(Number.NEGATIVE_INFINITY);
  t.set(hiSorted, size);
  for (let i = size - 1; i >= 1; i--) t[i] = Math.max(t[2 * i], t[2 * i + 1]);
  return t;
}

const stack = new Int32Array(64);
/**
 * The intervals i (sorted order) with `lo[i] ≤ q[1]` and a high end ≥ `q[0]`, ascending, into `out` (room for every
 * interval); their count.
 */
export function overlapping(lo: Float64Array, tree: Float64Array, q: Float64Array, out: Int32Array): number {
  let a = 0;
  let b = lo.length;
  while (a < b) {
    const m = (a + b) >> 1;
    if (lo[m] <= q[1]) a = m + 1;
    else b = m;
  }
  const last = a - 1;
  const size = tree.length >> 1;
  let n = 0;
  let sp = 0;
  if (last >= 0) stack[sp++] = 1;
  while (sp > 0) {
    const node = stack[--sp];
    IndexVisits.nodes++;
    if (tree[node] < q[0]) continue;
    const level = 31 - Math.clz32(node);
    if ((node - (1 << level)) * (size >> level) > last) continue;
    if (node >= size) {
      out[n++] = node - size;
      IndexVisits.intervals++;
      continue;
    }
    // The right child first, so the left is entered first: ascending.
    stack[sp++] = 2 * node + 1;
    stack[sp++] = 2 * node;
  }
  return n;
}
