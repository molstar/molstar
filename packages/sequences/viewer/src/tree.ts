/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tree beside a `rows` layer: x by branch length, y in display rows per layout (collapsed clades as wedge rows),
// the indices a band of rows draws from, the forgiving hit of a node, and the rows a node stands for.
// Headless; the rows renderer's prepared state (`TreeRows`) holds the geometry, the layout's `RowModel` the display
// rows. Its draw is the rows renderer's (RN3).

import type { SeqSet, Tree } from "@molstar/seq-core/model";
import type { AxisView } from "./axis.js";
import type { Hit } from "./hit.js";
import type { RowModel } from "./layout.js";

/** What the tree reads of the rows renderer's prepared state (its `RowsState` has these). */
export interface TreeRows {
  readonly set: SeqSet;
  readonly alignment: number;
  /** Lane offset → alignment row; null: a lazy run (`order`, else the offset itself). */
  readonly rowOf: Int32Array | null;
  readonly order: Int32Array | null;
  /** Labels of rows not in memory (a streamed alignment's), once served. */
  readonly source: { readonly remote: { readonly labels: ReadonlyMap<number, string> } | null } | null;
  readonly tree: TreeGeo | null;
}

/** A lane offset's alignment row. */
export const rowAt = (s: TreeRows, k: number) => (s.rowOf ? s.rowOf[k] : s.order ? s.order[k] : k);

/**
 * A tree beside the rows: x by branch length from the root (0..1); per layout, y in display rows (a leaf at its
 * row's center, an inner node between its children's extremes, a collapsed clade at its wedge's row), the display rows
 * each node stands for, and two indices for drawing a band of rows: nodes by y, and connectors by their top.
 */
export interface TreeGeo {
  readonly t: Tree;
  readonly x: Float32Array;
  readonly y: Float32Array;
  readonly lo: Float32Array;
  readonly hi: Float32Array;
  readonly first: Int32Array;
  readonly last: Int32Array;
  readonly state: Uint8Array;
  byY: Int32Array;
  ys: Float32Array;
  conn: Int32Array;
  maxHi: Float32Array;
  rowDisplay: Int32Array;
  /** Per collapsed clade (a wedge), its row's label: its first row's and how many more ("SYN001 +19"); made per layout. */
  wedgeText: Map<number, string>;
  /** The layout version laid out for; -1: none. */
  at: number;
}
/** `TreeGeo.state`: shown, hidden in a collapsed clade, a collapsed clade's wedge. */
export const TreeNode = { shown: 0, hidden: 1, wedge: 2 } as const;

export function treeOf(t: Tree): TreeGeo {
  const n = t.count;
  const x = new Float32Array(n);
  let depth = 0;
  for (let k = 0; k < n; k++) {
    const p = t.parent[k];
    x[k] = p < 0 ? 0 : x[p] + (Number.isNaN(t.branchLength[k]) ? 1 : t.branchLength[k]);
    depth = Math.max(depth, x[k]);
  }
  for (let k = 0; k < n; k++) x[k] /= depth || 1;
  const f = () => new Float32Array(n);
  const i = () => new Int32Array(n);
  return {
    t,
    x,
    y: f(),
    lo: f(),
    hi: f(),
    first: i(),
    last: i(),
    state: new Uint8Array(n),
    byY: i(),
    ys: f(),
    conn: i(),
    maxHi: f(),
    rowDisplay: new Int32Array(0),
    wedgeText: new Map(),
    at: -1,
  };
}

/**
 * Lays the tree out for the layout's rows, on a new version only. The work is apart: its sorts' closures would give
 * this check a context to allocate per call.
 */
export function layoutTree(s: TreeRows, r: RowModel) {
  if (s.tree && s.tree.at !== r.version) relayout(s, s.tree, r);
}
function relayout(s: TreeRows, g: TreeGeo, r: RowModel) {
  g.at = r.version;
  const { t, y, lo, hi, first, last, state } = g;
  const rows = s.set.alignments[s.alignment]?.rowCount ?? 0;
  if (g.rowDisplay.length !== rows) g.rowDisplay = new Int32Array(rows);
  const rd = g.rowDisplay.fill(-1);
  state.fill(TreeNode.shown);
  g.wedgeText.clear();
  lo.fill(Number.POSITIVE_INFINITY);
  hi.fill(Number.NEGATIVE_INFINITY);
  first.fill(0x7fffffff);
  last.fill(-1);
  for (let d = 0; d < r.count; d++) {
    const c = r.node ? r.node[d] : -1;
    if (c >= 0) {
      state[c] = TreeNode.wedge;
      y[c] = d + 0.5;
      first[c] = last[c] = d;
      state.fill(TreeNode.hidden, c + 1, t.cladeEnd[c]);
      let n = 0;
      let lead = -1;
      for (let k = c; k < t.cladeEnd[c]; k++)
        if (t.leafRow[k] >= 0) {
          if (lead < 0) lead = t.leafRow[k];
          n++;
        }
      g.wedgeText.set(c, n > 1 ? `${rowLabel(s, lead)} +${n - 1}` : rowLabel(s, lead));
      continue;
    }
    const k = (r.lane ? r.lane[d] : r.laneStart + d) - r.laneStart;
    const row = rowAt(s, k);
    if (row >= 0 && row < rows) rd[row] = d;
  }
  for (let k = t.count - 1; k >= 0; k--) {
    if (state[k] === TreeNode.hidden) continue;
    const row = t.leafRow[k];
    if (row >= 0) {
      const d = row < rows ? rd[row] : -1;
      y[k] = d >= 0 ? d + 0.5 : Number.NaN;
      if (d >= 0) first[k] = last[k] = d;
    } else if (state[k] === TreeNode.shown) y[k] = lo[k] <= hi[k] ? (lo[k] + hi[k]) / 2 : Number.NaN;
    const p = t.parent[k];
    if (p < 0 || Number.isNaN(y[k])) continue;
    lo[p] = Math.min(lo[p], y[k]);
    hi[p] = Math.max(hi[p], y[k]);
    first[p] = Math.min(first[p], first[k]);
    last[p] = Math.max(last[p], last[k]);
  }
  // Nodes by y; connectors (inner nodes with children shown) by their top, with the running maximum of bottoms.
  const nodes: number[] = [];
  const conns: number[] = [];
  for (let k = 0; k < t.count; k++) {
    if (state[k] === TreeNode.hidden || Number.isNaN(y[k])) continue;
    nodes.push(k);
    if (state[k] === TreeNode.shown && t.leafRow[k] < 0 && lo[k] < hi[k]) conns.push(k);
  }
  nodes.sort((a, b) => y[a] - y[b]);
  conns.sort((a, b) => lo[a] - lo[b]);
  g.byY = Int32Array.from(nodes);
  g.ys = Float32Array.from(nodes, (k) => y[k]);
  g.conn = Int32Array.from(conns);
  g.maxHi = new Float32Array(conns.length);
  for (let i = 0; i < conns.length; i++) g.maxHi[i] = Math.max(i > 0 ? g.maxHi[i - 1] : -1, hi[conns[i]]);
}

/** The tree's box in the header: x from its left over its width (header px), by its inset. */
export const TREE_INSET = 72;
export const treeX = (g: TreeGeo, k: number, left: number, width: number) => left + g.x[k] * width;
/** Nodes of `g.byY` from the first with y ≥ y0 (binary search). */
export function firstByY(g: TreeGeo, y0: number): number {
  let lo = 0;
  let hi = g.ys.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (g.ys[mid] < y0) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}
/** Connectors from the last whose top is below y1, back while the running bottom reaches y0: visit(k). */
export function lastConn(g: TreeGeo, y1: number): number {
  let lo = 0;
  let hi = g.conn.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (g.lo[g.conn[mid]] < y1) lo = mid + 1;
    else hi = mid;
  }
  return lo - 1;
}

/** The tree's hitbox: a node within 6 px; the header's hit reaches 8 px past the tree's inset. */
const TREE_HIT = 6;
const TREE_REACH = 8;

/**
 * The tree node under the point (`out.point`), a forgiving hitbox as the spike's: within the tree's inset + 8
 * px, the nearest branch (from its parent's x to its own, a wedge's triangle included) by its distance in y, or
 * connector by its distance in x, within 6 px. A leaf's branch is its row: nearest to a leaf, no node. `v.top` is display row 0's y and `v.height`
 * a row's (hit.ts sets them). Headless and allocation-free (the tree laid out for the layout).
 */
const nearest = new Float64Array(1);
export function hitTree(v: AxisView, s: TreeRows, out: Hit): boolean {
  const x = out.point[0];
  const y = out.point[1];
  const g = s.tree;
  const r = v.axis?.rows;
  const h = v.height;
  if (!g || !r || out.displayRow < 0 || !(h > 0) || x > TREE_INSET + TREE_REACH) return false;
  layoutTree(s, r);
  const w = TREE_INSET - TREE_LEFT - TREE_RIGHT;
  const top = v.top;
  // The point in display rows, and the reach in rows; positions and both searches inline (`firstByY`, `lastConn`): a
  // fractional double passed to or returned from a call V8 does not inline is boxed.
  const ry = (y - top) / h;
  const reach = TREE_HIT / h;
  const leaves = h >= 1;
  let best = -1;
  // The best distance in a typed array: a local phi of a constant and doubles may be kept boxed (16 bytes a node).
  nearest[0] = TREE_HIT;
  let first = 0;
  for (let hi = g.ys.length; first < hi; ) {
    const mid = (first + hi) >> 1;
    if (g.ys[mid] < ry - reach) first = mid + 1;
    else hi = mid;
  }
  for (let i = first; i < g.byY.length && g.ys[i] <= ry + reach; i++) {
    const k = g.byY[i];
    const p = g.t.parent[k];
    // Rows thinner than a pixel draw no leaves' branches (canvas/rows.ts): none to hit.
    if (p < 0 || (!leaves && g.t.leafRow[k] >= 0)) continue;
    const x0 = TREE_LEFT + g.x[p] * w - 2;
    const x1 = TREE_LEFT + g.x[k] * w + (g.state[k] === TreeNode.wedge ? 12 : 2);
    if (x < x0 || x > x1) continue;
    const d = Math.abs(y - (top + g.y[k] * h));
    if (d < nearest[0]) {
      best = k;
      nearest[0] = d;
    }
  }
  const slack = 2 / h;
  let last = 0;
  for (let hi = g.conn.length; last < hi; ) {
    const mid = (last + hi) >> 1;
    if (g.lo[g.conn[mid]] < ry + slack) last = mid + 1;
    else hi = mid;
  }
  for (let i = last - 1; i >= 0 && g.maxHi[i] >= ry - slack; i--) {
    const k = g.conn[i];
    if (g.hi[k] < ry - slack) continue;
    const d = Math.abs(x - (TREE_LEFT + g.x[k] * w));
    if (d < nearest[0]) {
      best = k;
      nearest[0] = d;
    }
  }
  if (best < 0 || g.t.leafRow[best] >= 0) return false;
  out.treeNode = best;
  return true;
}
export const TREE_LEFT = 6;
export const TREE_RIGHT = 14;

/** The display rows a tree node stands for, into out [first, last]; -1 when none shown. */
export function treeSpan(s: TreeRows, node: number, out: Int32Array) {
  const g = s.tree;
  out[0] = out[1] = -1;
  if (!g || node < 0 || node >= g.t.count || g.last[node] < 0) return;
  out[0] = g.first[node];
  out[1] = g.last[node];
}

/** A row's label: its alignment's (in memory), else its source's once served. */
export function rowLabel(s: TreeRows, row: number): string {
  const a = s.set.alignments[s.alignment];
  if (!a || row < 0) return "";
  const id = a.labelOf(row as never);
  if (id >= 0) return s.set.strings.get(id);
  return s.source?.remote?.labels.get(row) ?? "";
}
