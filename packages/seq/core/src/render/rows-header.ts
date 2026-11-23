/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { fontSize } from './config.js';
import type { LayoutContext } from './layer-type.js';
import { addName } from './row-names.js';
import { elide } from './rows-columns.js';
import { type Clades, columnsRoom, isGap } from './rows-model.js';
import { NoNumber, type RowsLay, type RowsPrep } from './rows-lay.js';

// The rows layer's header in the label column (_spec/README.md "Render", The rows layer: the header; The tree's look):
// the tree's segments, wedges and nodes from its kept geometry, the labels, the first and last residue numbers.

/** The tree's toggle glyph's side, CSS px (a narrow frame's collapsed tree). */
export const TogglePx = 12;

/** Whether a page point is on the tree's toggle (its glyph, a px around). */
export function inToggle(lay: RowsLay, x: number, y: number): boolean {
  return (
    lay.toggleX >= 0 &&
    x >= lay.toggleX - 1 &&
    x < lay.toggleX + TogglePx + 1 &&
    y >= lay.toggleY - 1 &&
    y < lay.toggleY + TogglePx + 1
  );
}

/**
 * The header in the label column: the tree (its own width, the label column widened by it), the labels (cut with an
 * ellipsis past their room), the first residue numbers before the data and the last after it (a right gutter).
 */
export function layoutHeader(lay: RowsLay, prep: RowsPrep, ctx: LayoutContext) {
  const { config, box, text } = ctx;
  const rs = lay.rowScale;
  const left = config.inset;
  const right = box.x - config.labelGap;
  const treeW =
    lay.clades && prep.order === 'tree' && !lay.figure
      ? Math.max(0, Math.min(ctx.treePx ?? prep.treeWidth, right - left - 24))
      : 0;
  // an opened tree keeps its toggle at the band's top with the tree beside it (pressed: it folds the tree back)
  const open = ctx.treeOpen === true && treeW > TogglePx + 4;
  const shift = open ? TogglePx + 4 : 0;
  lay.toggleOpen = open;
  lay.treeX0 = left + shift;
  lay.treeX1 = left + treeW;
  // the tree inset by half a stroke each side (`BranchPx`): a branch is centred on its line, and the header's clip is
  // the column's edge, so a vertical on the edge drew half its width
  if (treeW > BranchPx) layoutTree(lay, lay.clades!, left + shift + BranchPx / 2, treeW - shift - BranchPx);
  // a tree a narrow frame collapsed: its toggle in its place, at the band's top (a click or Shift+T reopens it)
  const toggle = ctx.treePx === 0 && lay.clades && prep.order === 'tree' && !lay.figure ? TogglePx + 4 : 0;
  if (toggle || open) {
    lay.toggleX = left;
    lay.toggleY = box.y + lay.pinnedHeight;
  }
  if (rs < config.rows.labelMinPx) return;
  // a row's label is data, as an item's: the items' font (a Pfam label, NAME/1-147, reads whole in the column)
  lay.labelFont = config.itemFont;
  lay.labelX = left + treeW + (treeW ? 4 : 0) + toggle;
  const numbers = prep.numbers === 'gutters';
  const numberW = numbers ? text.measure(config.itemFont, '00000') : 0;
  lay.numberX = right;
  lay.numberRightX = box.x + box.width + config.labelGap;
  // the left row columns between the label and the first numbers (the label column was widened by them)
  const leftRoom = columnsRoom(prep.columns, lay.cols.width, 'left', config.rows.columns.gapPx);
  lay.labelMax = Math.max(0, right - numberW - 4 - lay.labelX - leftRoom);
  const size = fontSize(config.itemFont);
  const add = (i: number) => {
    const e = lay.display.rowAt(i);
    const y = lay.rowY(i) + rs / 2 + size * 0.35;
    let t: string;
    if (e < 0) {
      // a summary's count is what it says (EC-100): the first label is cut, never ` +n`
      const m = lay.display.members.get(i);
      const more = lay.summaryMore(i);
      const head = m
        ? elide(text, lay.labelFont, lay.label(m[0]), lay.labelMax - text.measure(lay.labelFont, more))
        : '';
      t = m ? (head ? `${head}${more}` : elide(text, lay.labelFont, more.trim(), lay.labelMax)) : '';
    } else {
      t = lay.label(e);
      if (lay.display.pinned && i >= lay.display.pinned && config.rows.pinInPlace && isPinned(lay, e)) t = `· ${t}`;
      t = elide(text, lay.labelFont, t, lay.labelMax);
    }
    addName(text, lay.labelFont, t, y, lay.labelText, lay.labelY, lay.labelW);
    if (!numbers || e < 0 || !endNumbers(lay, e)) return;
    const n = String(lay.numFirst);
    lay.numberText.push(n);
    lay.numberY.push(y);
    lay.numberW.push(text.measure(config.itemFont, n));
    // the last residue number in the right gutter (a figure has none: its lines repeat the rows)
    if (!lay.figure) {
      const s = String(lay.numLast);
      const w = text.measure(config.itemFont, s);
      lay.lastText.push(s);
      lay.lastY.push(y);
      lay.lastW.push(w);
      if (w > lay.lastMaxW) lay.lastMaxW = w;
    }
  };
  for (let i = 0; i < lay.display.pinned; i++) add(i);
  for (let i = lay.first; i < lay.last; i++) add(i);
}

function isPinned(lay: RowsLay, row: number): boolean {
  for (let i = 0; i < lay.display.pinned; i++) if (lay.display.rowAt(i) === row) return true;
  return false;
}

/** The tree's stroke, px (`drawHeader`'s `branchWidth`). */
const BranchPx = 1;

/** A wedge's least length, px: a clade of zero branch lengths still shows one. */
const WedgeMinPx = 6;

/**
 * The least height a clade's rows take for its branch to be drawn, px (_spec/README.md "Render", The rows layer, The
 * tree's look): below it a clade's branches would sit under 2 px apart, and 10,000 leaves in a band read as a solid
 * slab; the tree is drawn down to the clades that resolve.
 */
const TreeMinPx = 4;

/**
 * The tree beside the rows from its kept geometry: only what crosses the band is kept (segments, each with the node a
 * hit on it toggles: an inner node's own vertical and the branch to it; a leaf's branch none), the wedges of folded
 * clades and the inner nodes a hit reads.
 */
function layoutTree(lay: RowsLay, clades: Clades, x0: number, w: number) {
  treeGeometry(lay, clades);
  const n = clades.nodes;
  const { coord, nodeX, inFold, folded, box } = lay;
  const rs = lay.rowScale;
  const P = lay.display.pinned;
  // a display row's y: the scrolled band's (the tree's rows are in the order, never the pinned band's copies)
  const y0 = box.y + lay.pinnedHeight - P * rs - lay.bandTop;
  const top = box.y + lay.pinnedHeight;
  const bottom = box.y + box.height;
  const inner = (k: number) => clades.hi[k] - clades.lo[k] > 1;
  // a clade drawn only when its rows take `TreeMinPx` (a node's segments, the branch to it, a folded clade's wedge)
  const least = TreeMinPx / rs;
  const span = lay.spanRows;
  for (let k = 0; k < n; k++) {
    if (inFold[k] || Number.isNaN(coord[k]) || span[k] < least) continue;
    const x = x0 + nodeX[k] * w;
    const yk = y0 + coord[k] * rs;
    if (folded[k]) {
      if (yk + rs / 2 >= top && yk - rs / 2 <= bottom) {
        // to its deepest leaf's x (a few px at least, so a clade of zero lengths still shows one); the hit's box the
        // row's middle 0.8, the drawn wedge's far side half the row (`drawWedges`)
        const far = Math.min(x0 + w, Math.max(x0 + lay.deepX[k] * w, x + WedgeMinPx));
        // a node within `WedgeMinPx` of the tree's edge: the wedge starts that far before it (its hit box too)
        lay.wedges.push(Math.min(x, far - WedgeMinPx), yk - rs * 0.4, far, yk + rs * 0.4, k);
        lay.nodes.push(k, x, yk);
      }
      continue;
    }
    const kids = clades.children[k];
    // the vertical spans the node and its drawn children only (none drawn: none; zoomed out, a child too small to
    // resolve has no branch to meet)
    let lo = yk;
    let hi = yk;
    let drawn = false;
    for (const c of kids) {
      if (Number.isNaN(coord[c]) || span[c] < least) continue;
      const y = y0 + coord[c] * rs;
      drawn = true;
      if (y < lo) lo = y;
      if (y > hi) hi = y;
      // the child's branch: from this node's x to the child's, at the child's y
      if (y >= top - 1 && y <= bottom + 1) {
        lay.tree.push(x, y, x0 + nodeX[c] * w, y);
        lay.treeNode.push(inner(c) ? c : -1);
        lay.treeClade.push(c);
      }
    }
    if (drawn && hi >= top && lo <= bottom) {
      lay.tree.push(x, Math.max(lo, top), x, Math.min(hi, bottom));
      lay.treeNode.push(inner(k) ? k : -1);
      lay.treeClade.push(k);
      if (yk >= top && yk <= bottom) lay.nodes.push(k, x, yk);
    }
  }
}

/** A clade's first shown row (its summary row's), -1 when none is shown. */
function firstShown(lay: RowsLay, clades: Clades, k: number): number {
  for (let l = clades.lo[k]; l < clades.hi[k]; l++) {
    const r = clades.leafRow[l];
    if (r >= 0 && lay.display.placeRow(r) >= 0) return r;
  }
  return -1;
}

/**
 * The residue numbers at the ends of a data row's part in view (_spec/README.md "Render", The rows layer: the header):
 * the numbers of its first and last residue in the columns in view, into `numFirst` and `numLast`; a gap is never
 * counted; false for a row with no residue in view. At the fit view they are the row's own (kept: counted once a data
 * or a source read), a pan or a zoom counts the residues before and in the columns in view.
 */
function endNumbers(lay: RowsLay, row: number): boolean {
  const held = lay.held;
  let codes: Uint8Array | undefined;
  let base = 0;
  let first: number;
  let C: number;
  if (held) {
    codes = held.codes;
    C = held.columns;
    base = row * C;
    first = held.first[row];
  } else {
    codes = lay.lazyCodes.get(row);
    if (!codes) return false;
    C = codes.length;
    first = lay.lazyFirst.get(row) ?? 1;
  }
  const c0 = lay.figure ? 0 : lay.viewC0;
  const c1 = lay.figure ? C : Math.min(lay.viewC1, C);
  if (c1 <= c0) return false;
  if (c0 === 0 && c1 === C) {
    const last = ownLast(lay, row, codes, base, C, first);
    if (last === NoNumber) return false;
    lay.numFirst = first;
    lay.numLast = last;
    return true;
  }
  let before = 0;
  for (let j = 0; j < c0; j++) if (!isGap(codes[base + j])) before++;
  let inView = 0;
  for (let j = c0; j < c1; j++) if (!isGap(codes[base + j])) inView++;
  if (!inView) return false;
  lay.numFirst = first + before;
  lay.numLast = first + before + inView - 1;
  return true;
}

/** A row's own last residue number: its first plus its residues less one (`NoNumber` for a row of gaps), kept. */
function ownLast(lay: RowsLay, row: number, codes: Uint8Array, base: number, C: number, first: number): number {
  const held = lay.held;
  if (held) {
    if (lay.lastOf !== held) {
      lay.lastOf = held;
      const out = new Int32Array(held.rows);
      for (let r = 0; r < held.rows; r++) {
        let n = 0;
        for (let j = 0; j < C; j++) if (!isGap(held.codes[r * C + j])) n++;
        out[r] = n ? held.first[r] + n - 1 : NoNumber;
      }
      lay.lasts = out;
    }
    return lay.lasts[row];
  }
  // a lazy row's: counted once when its tile is read, kept (bounded: let go past 4,096 rows)
  const kept = lay.lazyLasts.get(row);
  if (kept !== undefined) return kept;
  let n = 0;
  for (let j = 0; j < C; j++) if (!isGap(codes[base + j])) n++;
  const last = n ? first + n - 1 : NoNumber;
  if (lay.lazyLasts.size >= 4096) lay.lazyLasts.clear();
  lay.lazyLasts.set(row, last);
  return last;
}

/**
 * The tree's geometry in display rows, kept while the display is (_spec/README.md "Render", The rows layer: nothing
 * sized by the nodes made a pan): each node's x as a share of the tree's width (branch lengths; depth when none), its
 * y as a display row (a leaf at its row's place in the order, a pinned one in place; a folded clade at its summary
 * row; an inner node between its first and last child), whether it is folded or inside a fold.
 */
function treeGeometry(lay: RowsLay, clades: Clades) {
  const d = lay.display;
  if (lay.treeOf === d.version && lay.treeClades === clades) return;
  lay.treeOf = d.version;
  lay.treeClades = clades;
  const n = clades.nodes;
  if (lay.coord.length < n) {
    lay.coord = new Float64Array(n);
    lay.nodeX = new Float64Array(n);
    lay.deepX = new Float64Array(n);
    lay.spanRows = new Float64Array(n);
    lay.inFold = new Uint8Array(n);
    lay.folded = new Uint8Array(n);
  }
  const { coord, nodeX, deepX, inFold, folded } = lay;
  let maxD = 0;
  for (let k = 0; k < n; k++) if (clades.distance[k] > maxD) maxD = clades.distance[k];
  let maxDepth = 0;
  for (let k = 0; k < n; k++) if (clades.depth[k] > maxDepth) maxDepth = clades.depth[k];
  for (let k = 0; k < n; k++) nodeX[k] = maxD > 0 ? clades.distance[k] / maxD : clades.depth[k] / Math.max(1, maxDepth);
  // a node's children come after it: its deepest leaf's x from theirs
  for (let k = n - 1; k >= 0; k--) {
    const kids = clades.children[k];
    let d = nodeX[k];
    for (const c of kids) if (deepX[c] > d) d = deepX[c];
    deepX[k] = d;
  }
  // folded: a clade shown as a summary row (one folded to its one shown row is that row)
  folded.fill(0, 0, n);
  for (const i of d.members.keys()) {
    const e = d.rowAt(i);
    if (e < 0) folded[-1 - e] = 1;
  }
  for (let k = 0; k < n; k++) {
    const p = clades.parent[k];
    inFold[k] = p >= 0 && (inFold[p] || folded[p]) ? 1 : 0;
  }
  const { spanRows } = lay;
  for (let k = n - 1; k >= 0; k--) {
    coord[k] = Number.NaN;
    spanRows[k] = 0;
    if (inFold[k]) continue;
    const kids = clades.children[k];
    if (!kids.length || folded[k]) {
      const leafRow = kids.length ? -1 : clades.leafRow[clades.lo[k]];
      const row = leafRow >= 0 ? leafRow : firstShown(lay, clades, k);
      const at = row >= 0 ? d.placeRow(row) : -1;
      if (at >= 0) {
        coord[k] = at + 0.5;
        spanRows[k] = 1;
      }
      continue;
    }
    let lo = Number.NaN;
    let hi = Number.NaN;
    for (const c of kids) {
      const y = coord[c];
      if (Number.isNaN(y)) continue;
      if (Number.isNaN(lo) || y < lo) lo = y;
      if (Number.isNaN(hi) || y > hi) hi = y;
      spanRows[k] += spanRows[c];
    }
    if (!Number.isNaN(lo)) coord[k] = (lo + hi) / 2;
  }
}
