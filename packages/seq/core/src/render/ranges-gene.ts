/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { RangesList } from '@molstar/seqx-core/data';
import type { AxisTransform } from './axis.js';
import { grown } from './draw-lists.js';
import { type LayerSource, toSection } from './layer-data.js';
import { mappedPieces, type RangesLay } from './ranges-lay.js';

// The ranges layer's gene line (`line: 'gene'`; _spec/README.md "Render", Layer types, The gene line): an item's rows
// are its pieces (a transcript's exons, a BED12's blocks) in one lane, laid by the item's whole extent (`itemStart`,
// `itemEnd`, else the hull of its rows in view); each piece a box, full height between `thick_start` and `thick_end`
// (the CDS) and thin outside (the UTRs), all full when the item has no thick part (NaN) and all thin when it is empty
// (`thick_start = thick_end`: a non-coding BED12); a hairline across the item's extent under the boxes, with chevrons
// by its strand (`cdir`: 1 pointing right, 2 left; the connectors' arrays: a hit on it is the item's, the hover outlines it). A span row (no piece in its
// tile) is no box, only the item. Below a position a pixel too (the boxes are never merged into runs).

/** Where an item's lane is laid: the first row of its key (`repOf`) holds the item's extent (`gStart`, `gEnd`). */
export interface GeneScratch {
  repOf: Int32Array;
  gStart: Float64Array;
  gEnd: Float64Array;
  /** The rows by key, cleared each layout. */
  readonly groups: Map<string | number, number>;
  /** The last table's columns (rows of one table share its fields). */
  fields: RangesList['fields'] | undefined;
  ts: Column | undefined;
  te: Column | undefined;
  strand: Column | undefined;
}

export function newGeneScratch(): GeneScratch {
  return {
    repOf: new Int32Array(64),
    gStart: new Float64Array(64),
    gEnd: new Float64Array(64),
    groups: new Map(),
    fields: undefined,
    ts: undefined,
    te: undefined,
    strand: undefined,
  };
}

/**
 * Groups the `n` rows in view by item: each row's first row of its key (`repOf`), the item's extent on that row, and
 * the first rows into `order` (returned: their count), so the lanes are packed an item at a time.
 */
export function geneGroups(lay: RangesLay, n: number, order: Int32Array): number {
  const g = lay.gene;
  if (g.repOf.length < n) {
    g.repOf = grown(g.repOf, n);
    g.gStart = grown(g.gStart, n);
    g.gEnd = grown(g.gEnd, n);
  }
  g.groups.clear();
  const { keys, itemStarts, itemEnds } = lay.rows;
  let count = 0;
  for (let j = 0; j < n; j++) {
    const rep = g.groups.get(keys[j]);
    const known = itemStarts[j] >= 0 && itemEnds[j] >= itemStarts[j];
    if (rep === undefined) {
      g.groups.set(keys[j], j);
      g.repOf[j] = j;
      g.gStart[j] = known ? itemStarts[j] : lay.starts[j];
      g.gEnd[j] = known ? itemEnds[j] : lay.ends[j];
      order[count++] = j;
      continue;
    }
    g.repOf[j] = rep;
    if (!known) {
      g.gStart[rep] = Math.min(g.gStart[rep], lay.starts[j]);
      g.gEnd[rep] = Math.max(g.gEnd[rep], lay.ends[j]);
    }
  }
  return count;
}

function numberAt(c: Column | undefined, i: number): number {
  if (!c || c instanceof CategoryColumn) return Number.NaN;
  return (c as ArrayLike<number>)[i];
}

/** Binds the scratch to row `j`'s table (its thick part's and strand's columns). */
function bind(lay: RangesLay, j: number) {
  const g = lay.gene;
  const fields = lay.rows.fields[j];
  if (fields === g.fields) return;
  g.fields = fields;
  g.ts = g.te = g.strand = undefined;
  if (!fields) return;
  for (let k = 0; k < fields.length; k++) {
    const name = fields[k].name;
    if (name === 'thick_start') g.ts = fields[k].column;
    else if (name === 'thick_end') g.te = fields[k].column;
    else if (name === 'strand') g.strand = fields[k].column;
  }
}

/**
 * The marks of the gene line: per row in a lane, its boxes (thick and thin, `mthin`), and per item its hairline with its
 * strand's direction on the screen (`cdir`). `next` makes a mark (the layer's own, its arrays grown).
 */
export function geneMarks(
  lay: RangesLay,
  axis: AxisTransform,
  src: LayerSource,
  n: number,
  left: number,
  right: number,
  next: (lay: RangesLay) => number,
): void {
  const g = lay.gene;
  const rows = lay.rows;
  const px = 1 / axis.dpr;
  for (let i = 0; i < n; i++) {
    const lane = lay.lanes[i];
    if (lane < 0) continue;
    bind(lay, i);
    const r = rows.row[i];
    if (g.repOf[i] === i) {
      const st = numberAt(g.strand, r);
      hairline(lay, axis, src, g.gStart[i], g.gEnd[i], lane, i, st > 0 ? 1 : st < 0 ? -1 : 0, left, right);
    }
    const s = lay.starts[i];
    const e = lay.ends[i];
    // a span row: the item, no piece here
    if (e <= s && rows.itemStarts[i] >= 0 && rows.itemStarts[i] < rows.itemEnds[i]) continue;
    const ts = numberAt(g.ts, r);
    const te = numberAt(g.te, r);
    if (Number.isNaN(ts) || Number.isNaN(te)) {
      boxes(lay, axis, src, s, e, lane, i, 0, px, next);
      continue;
    }
    if (te <= ts) {
      boxes(lay, axis, src, s, e, lane, i, 1, px, next);
      continue;
    }
    if (s < ts) boxes(lay, axis, src, s, Math.min(e, ts), lane, i, 1, px, next);
    if (Math.min(e, te) > Math.max(s, ts))
      boxes(lay, axis, src, Math.max(s, ts), Math.min(e, te), lane, i, 0, px, next);
    if (e > te) boxes(lay, axis, src, Math.max(s, te), e, lane, i, 1, px, next);
  }
}

/** The boxes of `[s, e)` of row `i` on the section, square, thin or not. */
function boxes(
  lay: RangesLay,
  axis: AxisTransform,
  src: LayerSource,
  s: number,
  e: number,
  lane: number,
  i: number,
  thin: number,
  px: number,
  next: (lay: RangesLay) => number,
) {
  const pieces = toSection(src, s, e, mappedPieces);
  for (let k = 0; k + 1 < pieces.length; k += 2) {
    axis.piecesIn(pieces[k], pieces[k + 1]);
    const ka = axis.ka;
    const kb = axis.kb;
    for (let a = ka; a <= kb; a++) {
      if (!axis.spanIn(a, pieces[k], pieces[k + 1])) continue;
      const x0 = axis.snap(axis.sx0);
      const j = next(lay);
      lay.mx0[j] = x0;
      lay.mx1[j] = Math.max(axis.snap(axis.sx1), x0 + px);
      lay.mlane[j] = lane;
      lay.mitem[j] = i;
      lay.mround[j] = 0;
      lay.mthin[j] = thin;
    }
  }
}

/** The hairline of an item's extent `[s, e)` (its data system's) in `lane`, clipped to `[left, right]`. */
function hairline(
  lay: RangesLay,
  axis: AxisTransform,
  src: LayerSource,
  s: number,
  e: number,
  lane: number,
  item: number,
  strand: number,
  left: number,
  right: number,
) {
  if (lay.calpha <= 0 || !(e > s)) return;
  const pieces = toSection(src, s, e, mappedPieces);
  if (!pieces.length) return;
  let lo = pieces[0];
  let hi = pieces[1];
  for (let k = 2; k + 1 < pieces.length; k += 2) {
    lo = Math.min(lo, pieces[k]);
    hi = Math.max(hi, pieces[k + 1]);
  }
  const xa = axis.x(lo);
  const xb = axis.x(hi);
  const x0 = Math.max(Math.min(xa, xb), left);
  const x1 = Math.min(Math.max(xa, xb), right);
  if (!(x1 > x0)) return;
  const j = lay.cn++;
  if (j >= lay.cx0.length) {
    lay.cx0 = grown(lay.cx0, j + 1);
    lay.cx1 = grown(lay.cx1, j + 1);
    lay.clane = grown(lay.clane, j + 1);
    lay.citem = grown(lay.citem, j + 1);
    lay.cdir = grown(lay.cdir, j + 1);
    lay.canchor = grown(lay.canchor, j + 1);
  }
  lay.canchor[j] = Math.min(xa, xb);
  lay.cx0[j] = x0;
  lay.cx1[j] = x1;
  lay.clane[j] = lane;
  lay.citem[j] = item;
  // on the screen: the strand read along the axis (a reversed axis turns it)
  const d = strand * axis.dir;
  lay.cdir[j] = d > 0 ? 1 : d < 0 ? 2 : 0;
}
