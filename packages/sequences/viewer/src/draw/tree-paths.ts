/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// (its `drawTree`, `drawTreeHover` and `treePaths`)
// The tree beside a rows layer, part of the renderer API: its branches, connectors and collapsed clades' wedges over
// display rows [first, last] (`v.top`: display row 0's y, `v.height` a row's, `v.left` the header's), as one stroked
// subpath and a filled triangle per wedge, in device px. Lines are a whole number of device px wide, centered on device
// px, so a copied band equals a redraw. Batches owned and handed over once.

import { firstByY, lastConn, layoutTree, TREE_INSET, TREE_LEFT, TREE_RIGHT, TreeNode, type TreeRows } from '../tree.js';
import { PathBatch, PathMode, PathVerb, SkinToken } from './backend.js';
import type { DrawView } from './view.js';

const paths = new PathBatch(32);
/** A segment's ends (x0, y0, x1, y1), device px: a double passed to a call is boxed. */
const seg = new Float64Array(4);

/** The segment in `seg` as a move and a line of the current subpath. */
function line(): void {
  let o = paths.verb(PathVerb.move);
  paths.coords[o] = seg[0];
  paths.coords[o + 1] = seg[1];
  o = paths.verb(PathVerb.line);
  paths.coords[o] = seg[2];
  paths.coords[o + 1] = seg[3];
}

/**
 * Nodes [n0, n1) (a clade: preorder) of the tree, stroked `wide` times the least line (one device px) in `ink`, `wide`
 * 1 centered on the edge's pixels, 2 on the pixel edge (the hover's, so it covers the tree's lines).
 */
function treePaths(
  v: DrawView,
  s: TreeRows,
  first: number,
  last: number,
  n0: number,
  n1: number,
  wide: number,
  ink: number,
) {
  const g = s.tree;
  const r = v.axis?.rows;
  if (!g || !r) return;
  layoutTree(s, r);
  const d = v.dpr;
  const lw = Math.max(1, Math.round(d)) | 0;
  const half = wide === 1 ? lw / 2 : 0;
  const h = v.height;
  const left = v.left + TREE_LEFT;
  const w = TREE_INSET - TREE_LEFT - TREE_RIGHT;
  paths.reset();
  paths.begin(ink, 255, PathMode.stroke, 2 * wide * lw);
  for (let i = firstByY(g, first); i < g.byY.length && g.ys[i] < last + 1; i++) {
    const k = g.byY[i];
    const p = g.t.parent[k];
    // Rows thinner than a device pixel: leaves' branches would fill the box; inner nodes and connectors read.
    if (k < n0 || k >= n1 || p < 0 || (h * d < 1 && g.t.leafRow[k] >= 0)) continue;
    seg[1] = seg[3] = Math.round((v.top + g.y[k] * h) * d) + half;
    seg[0] = Math.round((left + g.x[p] * w) * d) + half;
    seg[2] = Math.round((left + g.x[k] * w) * d) + half;
    line();
  }
  for (let i = lastConn(g, last + 1); i >= 0 && g.maxHi[i] > first; i--) {
    const k = g.conn[i];
    if (k < n0 || k >= n1 || g.hi[k] <= first) continue;
    seg[0] = seg[2] = Math.round((left + g.x[k] * w) * d) + half;
    // The connector's rows clamped to the band, in `seg` (a conditional over doubles keeps a boxed number). An end cut
    // at the band's edge runs a line width past it, for the band's clip to cut: ended on the edge, the line's half-pixel
    // centering left the edge's pixel row half covered in each band, a light gap where bands meet.
    seg[1] = g.lo[k];
    seg[3] = g.hi[k];
    const cutTop = seg[1] < first;
    const cutBottom = seg[3] > last + 1;
    if (cutTop) seg[1] = first;
    if (cutBottom) seg[3] = last + 1;
    seg[1] = Math.round((v.top + seg[1] * h) * d) + half - (cutTop ? lw : 0);
    seg[3] = Math.round((v.top + seg[3] * h) * d) + half + (cutBottom ? lw : 0);
    line();
  }
  const hh = (h < 15 ? h * 0.4 : 6) * d;
  for (let i = firstByY(g, first); hh >= d && i < g.byY.length && g.ys[i] < last + 1; i++) {
    const k = g.byY[i];
    if (k < n0 || k >= n1 || g.state[k] !== TreeNode.wedge) continue;
    const x = Math.round((left + g.x[k] * w) * d);
    const y = Math.round((v.top + g.y[k] * h) * d);
    paths.begin(ink, 255, PathMode.fill);
    let o = paths.verb(PathVerb.move);
    paths.coords[o] = x;
    paths.coords[o + 1] = y;
    o = paths.verb(PathVerb.line);
    paths.coords[o] = x + 10 * d;
    paths.coords[o + 1] = y - hh;
    o = paths.verb(PathVerb.line);
    paths.coords[o] = x + 10 * d;
    paths.coords[o + 1] = y + hh;
    paths.close();
  }
  v.backend.paths(paths);
}

/** The tree over display rows [first, last]: branches into the nodes on those rows, the connectors crossing them, wedges. */
export function drawTree(v: DrawView, s: TreeRows, first: number, last: number): void {
  if (s.tree) treePaths(v, s, first, last, 0, s.tree.t.count, 1, SkinToken.inkQuiet);
}
/** A hovered clade (`node` and the nodes under it) in the accent, twice the tree's width and centered on its lines. */
export function drawTreeHover(v: DrawView, s: TreeRows, node: number, first: number, last: number): void {
  const g = s.tree;
  if (g && node >= 0 && node < g.t.count) treePaths(v, s, first, last, node, g.t.cladeEnd[node], 2, SkinToken.accent);
}
