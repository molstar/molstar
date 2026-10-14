/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Guide trees (seqspec:model#alignments-and-trees, seqspec:model#providers): a distance matrix over an alignment's
// rows, and the two textbook builders over it, neighbor joining (Saitou and Nei 1987) and UPGMA (Sokal and Michener
// 1958). Pure, over typed arrays, generic over any alphabet (a column is a letter or none). Both are O(rows^3) in time
// and O(rows^2) in memory, so they are bounded (`MAX_TREE_ROWS`) and take an abort signal, which they honor between
// slices of work: a long build yields to the event loop every `SLICE` basic steps (`Pacer`, pace.ts) and throws the
// signal's reason when it was aborted.

import { Pacer } from './pace.js';

/** The most rows a guide tree is built over: the matrix takes rows^2 doubles (32 MB here) and the build rows^3 / 6 steps. */
export const MAX_TREE_ROWS = 2048;

/** Thrown before anything is allocated when a tree would take more than `MAX_TREE_ROWS` rows. */
export class TreeSizeError extends Error {
  readonly rows: number;
  constructor(rows: number) {
    super(`a guide tree over ${rows} rows takes ${rows}^2 distances; the limit is ${MAX_TREE_ROWS} rows.`);
    this.name = 'TreeSizeError';
    this.rows = rows;
  }
}

/**
 * A rooted tree as parallel arrays: nodes `0 .. leaves - 1` are the leaves, node k the input's row k; the rest are inner
 * nodes in creation order, the root last. `branchLength` is the edge to the parent (NaN at the root).
 */
export interface GuideTree {
  readonly leaves: number;
  readonly parent: Int32Array;
  readonly branchLength: Float64Array;
}

/**
 * The distance matrix of `rowCount` rows of `columnCount` columns in `codes` (row-major, 0 where a row has no letter, any
 * other value a letter): per pair, the share of the columns where both rows have a letter at which the letters differ.
 * A pair with no such column is 1 (the largest). Symmetric, zero on the diagonal, row-major `rowCount^2`.
 * `p` and `identity` are this one number.
 */
export async function alignmentDistances(
  codes: Uint8Array,
  rowCount: number,
  columnCount: number,
  signal?: AbortSignal,
): Promise<Float64Array> {
  if (rowCount > MAX_TREE_ROWS) throw new TreeSizeError(rowCount);
  const pace = new Pacer(signal);
  const d = new Float64Array(rowCount * rowCount);
  for (let i = 0; i < rowCount; i++) {
    const oi = i * columnCount;
    for (let j = i + 1; j < rowCount; j++) {
      const oj = j * columnCount;
      let both = 0;
      let same = 0;
      for (let c = 0; c < columnCount; c++) {
        const a = codes[oi + c];
        const b = codes[oj + c];
        if (a !== 0 && b !== 0) {
          both++;
          if (a === b) same++;
        }
      }
      const v = both === 0 ? 1 : (both - same) / both;
      d[i * rowCount + j] = v;
      d[j * rowCount + i] = v;
    }
    if (pace.due((rowCount - i) * columnCount)) await pace.yieldNow();
  }
  return d;
}

/** Compacts `slots[0 .. m)` without the entry at index `at`. */
function drop(slots: Int32Array, m: number, at: number) {
  slots.copyWithin(at, at + 1, m);
}

/**
 * Neighbor joining over `n` rows' distances (`dist`, row-major `n^2`; not changed): the pair with the smallest
 * Q(i, j) = (m - 2) d(i, j) - r_i - r_j is joined each step; ties go to the first pair in row order (a joined pair's
 * node takes the lower slot). The tree is the textbook's unrooted one with its last three nodes joined at the root, so
 * branch lengths are the textbook's, except that a negative one is 0. Rows 1 and 2 are a root with one and two leaves.
 */
export async function neighborJoining(dist: Float64Array, n: number, signal?: AbortSignal): Promise<GuideTree> {
  if (n > MAX_TREE_ROWS) throw new TreeSizeError(n);
  const pace = new Pacer(signal);
  const total = Math.max(1, n === 1 ? 1 : n < 3 ? n + 1 : 2 * n - 2);
  const parent = new Int32Array(total).fill(-1);
  const length = new Float64Array(total).fill(Number.NaN);
  if (n === 2) {
    parent[0] = parent[1] = 2;
    length[0] = length[1] = Math.max(0, dist[1]) / 2;
    return { leaves: 2, parent, branchLength: length };
  }
  if (n < 2) return { leaves: n, parent, branchLength: length };
  const d = Float64Array.from(dist);
  const slots = new Int32Array(n); // the matrix rows still in play, in row order
  const node = new Int32Array(n); // a slot's tree node
  const r = new Float64Array(n); // a slot's sum of distances to the others in play
  for (let i = 0; i < n; i++) {
    slots[i] = node[i] = i;
    for (let j = 0; j < n; j++) r[i] += d[i * n + j];
  }
  let m = n;
  let next = n;
  while (m > 3) {
    let best = Number.POSITIVE_INFINITY;
    let bi = 0;
    let bj = 1;
    for (let a = 0; a < m - 1; a++) {
      const i = slots[a];
      const row = i * n;
      const base = r[i];
      for (let b = a + 1; b < m; b++) {
        const j = slots[b];
        const q = (m - 2) * d[row + j] - base - r[j];
        if (q < best) {
          best = q;
          bi = a;
          bj = b;
        }
      }
    }
    const f = slots[bi];
    const g = slots[bj];
    const dfg = d[f * n + g];
    const bf = dfg / 2 + (r[f] - r[g]) / (2 * (m - 2));
    parent[node[f]] = parent[node[g]] = next;
    length[node[f]] = Math.max(0, bf);
    length[node[g]] = Math.max(0, dfg - bf);
    // The new node takes f's slot; g leaves.
    let sum = 0;
    for (let a = 0; a < m; a++) {
      const k = slots[a];
      if (k === f || k === g) continue;
      const v = (d[f * n + k] + d[g * n + k] - dfg) / 2;
      r[k] += v - d[f * n + k] - d[g * n + k];
      d[f * n + k] = d[k * n + f] = v;
      sum += v;
    }
    r[f] = sum;
    node[f] = next++;
    drop(slots, m, bj);
    m--;
    if (pace.due(m * m)) await pace.yieldNow();
  }
  const [x, y, z] = [slots[0], slots[1], slots[2]];
  const [dxy, dxz, dyz] = [d[x * n + y], d[x * n + z], d[y * n + z]];
  const root = next;
  parent[node[x]] = parent[node[y]] = parent[node[z]] = root;
  length[node[x]] = Math.max(0, (dxy + dxz - dyz) / 2);
  length[node[y]] = Math.max(0, (dxy + dyz - dxz) / 2);
  length[node[z]] = Math.max(0, (dxz + dyz - dxy) / 2);
  return { leaves: n, parent, branchLength: length };
}

/**
 * The tree rerooted at the midpoint of its longest leaf-to-leaf path (found by two sweeps, the farthest leaf of the farthest leaf from leaf 0): a new root node
 * (last) splits that path's middle edge, and every other node keeps its index, so leaf k is still row k. For the
 * unrooted tree `neighborJoining` gives, whose root is where the last three nodes happened to meet. Fewer than three
 * leaves, or leaves all at distance 0, are returned as they are.
 */
export function midpointRoot(t: GuideTree): GuideTree {
  const n = t.parent.length;
  if (t.leaves < 3) return t;
  // Undirected adjacency (compressed rows) with each edge's length.
  const start = new Int32Array(n + 1);
  for (let v = 0; v < n; v++)
    if (t.parent[v] >= 0) {
      start[v + 1]++;
      start[t.parent[v] + 1]++;
    }
  for (let v = 0; v < n; v++) start[v + 1] += start[v];
  const to = new Int32Array(start[n]);
  const fill = start.slice(0, n);
  for (let v = 0; v < n; v++) {
    const p = t.parent[v];
    if (p >= 0) {
      to[fill[v]++] = p;
      to[fill[p]++] = v;
    }
  }
  const edge = (a: number, b: number) => (t.parent[a] === b ? t.branchLength[a] : t.branchLength[b]);
  const dist = new Float64Array(n);
  const prev = new Int32Array(n);
  const stack = new Int32Array(n);
  /** Distances from `from` to every node into `dist`, and the node before each on its path into `prev`; the farthest leaf. */
  const sweep = (from: number): number => {
    let sp = 0;
    stack[sp++] = from;
    dist[from] = 0;
    prev[from] = -1;
    let far = from;
    while (sp > 0) {
      const v = stack[--sp];
      if (v < t.leaves && dist[v] > dist[far]) far = v;
      for (let k = start[v]; k < start[v + 1]; k++) {
        const w = to[k];
        if (w === prev[v]) continue;
        prev[w] = v;
        dist[w] = dist[v] + edge(v, w);
        stack[sp++] = w;
      }
    }
    return far;
  };
  const u = sweep(0);
  const v = sweep(u);
  // Every leaf at distance 0 (a homo-oligomer's identical chains): no path has a middle, and any root is as good.
  if (dist[v] <= 0) return t;
  const half = dist[v] / 2;
  // The path from v back to u: the first edge (a, b) whose far end `a` is at or past the middle.
  let a = v;
  while (dist[prev[a]] >= half && prev[a] !== u) a = prev[a];
  const b = prev[a];
  const w = edge(a, b);
  const x = Math.max(0, Math.min(w, half - dist[b])); // from b to the new root
  // Orient everything away from the new root `n`, which has b and a as children.
  const parent = new Int32Array(n + 1);
  const length = new Float64Array(n + 1);
  parent[n] = -1;
  length[n] = Number.NaN;
  for (const [c, l] of [
    [b, x],
    [a, w - x],
  ] as const) {
    parent[c] = n;
    length[c] = l;
  }
  let sp = 0;
  stack[sp++] = b;
  stack[sp++] = a;
  while (sp > 0) {
    const p = stack[--sp];
    for (let k = start[p]; k < start[p + 1]; k++) {
      const c = to[k];
      if (c === parent[p] || (p === a && c === b) || (p === b && c === a)) continue;
      parent[c] = p;
      length[c] = edge(p, c);
      stack[sp++] = c;
    }
  }
  return { leaves: t.leaves, parent, branchLength: length };
}

/**
 * UPGMA over `n` rows' distances (`dist`, row-major `n^2`; not changed): the closest pair of clusters is joined each
 * step (ties: the first pair in row order, the joined cluster taking the lower slot), at height half their distance, and
 * a cluster's distance to another is the average over their rows. An ultrametric binary tree; branch lengths are height
 * differences. One row is a root alone; two are a root with two leaves.
 */
export async function upgma(dist: Float64Array, n: number, signal?: AbortSignal): Promise<GuideTree> {
  if (n > MAX_TREE_ROWS) throw new TreeSizeError(n);
  const pace = new Pacer(signal);
  const total = Math.max(1, 2 * n - 1);
  const parent = new Int32Array(total).fill(-1);
  const length = new Float64Array(total).fill(Number.NaN);
  if (n < 2) return { leaves: n, parent, branchLength: length };
  const d = Float64Array.from(dist);
  const slots = new Int32Array(n);
  const node = new Int32Array(n);
  const size = new Float64Array(n).fill(1);
  const height = new Float64Array(n);
  for (let i = 0; i < n; i++) slots[i] = node[i] = i;
  let m = n;
  let next = n;
  while (m > 1) {
    let best = Number.POSITIVE_INFINITY;
    let bi = 0;
    let bj = 1;
    for (let a = 0; a < m - 1; a++) {
      const row = slots[a] * n;
      for (let b = a + 1; b < m; b++) {
        const v = d[row + slots[b]];
        if (v < best) {
          best = v;
          bi = a;
          bj = b;
        }
      }
    }
    const f = slots[bi];
    const g = slots[bj];
    const h = best / 2;
    parent[node[f]] = parent[node[g]] = next;
    length[node[f]] = Math.max(0, h - height[f]);
    length[node[g]] = Math.max(0, h - height[g]);
    const [sf, sg] = [size[f], size[g]];
    for (let a = 0; a < m; a++) {
      const k = slots[a];
      if (k === f || k === g) continue;
      d[f * n + k] = d[k * n + f] = (sf * d[f * n + k] + sg * d[g * n + k]) / (sf + sg);
    }
    size[f] = sf + sg;
    height[f] = h;
    node[f] = next++;
    drop(slots, m, bj);
    m--;
    if (pace.due(m * m)) await pace.yieldNow();
  }
  return { leaves: n, parent, branchLength: length };
}
