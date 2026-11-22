/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArrayRange } from '@molstar/core/data/util/sort';
import { grown } from './draw-lists.js';
import { itemArc, itemHit } from './items.js';
import { type LayerSource, toSection } from './layer-data.js';
import type { LabelSink, LayoutContext } from './layer-type.js';
import { LabelLane } from './sizing.js';
import { labelsInside } from './range-labels.js';
import type { RangesLay } from './ranges-lay.js';

// A ranges layer's labels beside its marks (_spec/README.md "Render", Labels): placed by priority in the track's label
// lane, a run labeled once, a leader where a label moved; or every shown item's label to a layout's sink (rings).

/**
 * Labels inside if they fit, else beside (_spec/README.md "Render", Labels; range-labels.ts): in the label lane above
 * the track, at its item's start, a leader when it had to move (`LabelGapPx` between two labels there).
 */
const LabelGapPx = 6;

/**
 * Every shown item's label to a layout's sink (_spec/README.md "Render", Layer types: `labels`), an item once: its
 * extent on the section's system (its image's first and last edge: past each other over a cyclic origin), its element
 * and its hit, ordered by its band's place (the layer's order).
 */
export function sinkLabels(lay: RangesLay, src: LayerSource, sink: LabelSink, order: number) {
  // an item's rows (its pieces) by key, on the section's system
  const byKey = new Map<string | number, { text: string; row: number; pieces: number[] }>();
  for (let i = 0; i < lay.items; i++) {
    const text = lay.labels[i];
    if (!text || lay.lanes[i] < 0) continue;
    const key = lay.keys[i];
    let g = byKey.get(key);
    if (!g) byKey.set(key, (g = { text, row: i, pieces: [] }));
    if (src.path) g.pieces.push(...toSection(src, lay.starts[i], lay.ends[i], []));
    else g.pieces.push(lay.starts[i], lay.ends[i]);
  }
  const base = { layer: lay.key, section: lay.section, system: lay.system };
  const length = lay.axis.cyclic ? lay.axis.length : undefined;
  for (const [key, g] of byKey) {
    const [s, e] = itemArc(g.pieces, length);
    const hit = itemHit(base, lay, key, lay.starts[g.row], lay.ends[g.row], g.row);
    sink.label(g.text, s, e, { owner: lay.owner, key }, hit, order);
  }
}

/** The candidates' order (the longer item first, then the earlier): a comparator of the module, reading these. */
let SortStarts: readonly number[] = [];
let SortEnds: readonly number[] = [];
function longerFirst(o: ArrayLike<number>, a: number, b: number): number {
  const x = o[a];
  const y = o[b];
  return SortEnds[y] - SortStarts[y] - (SortEnds[x] - SortStarts[x]) || x - y;
}

/** The label lane of a layer laid out alone (no track's: a test's, an export's own call). */
const OwnLane = new LabelLane();

/**
 * Each label inside the first of its marks it fits whole; an item whose label fits none of its marks in view goes beside,
 * in the track's label lane at its leftmost mark's x (moved right past a label placed before it, the track's earlier
 * layers' first, a leader drawn from its own lane), placed by one priority (the layer's order, then the longer item);
 * one the lane cannot hold goes inside cut, else is dropped. Allocation-free once warm: per layout, never on a hover.
 */
export function placeLabels(lay: RangesLay, ctx: LayoutContext) {
  const n = lay.items;
  lay.left = grown(lay.left, n);
  const left = lay.left;
  left.fill(Number.POSITIVE_INFINITY, 0, n);
  lay.right = grown(lay.right, n);
  lay.right.fill(Number.NEGATIVE_INFINITY, 0, n);
  for (let j = 0; j < lay.hn; j++) {
    const i = lay.hitem[j];
    if (i < 0) continue;
    if (lay.hx0[j] < left[i]) left[i] = lay.hx0[j];
    if (lay.hx1[j] > lay.right[i]) lay.right[i] = lay.hx1[j];
  }
  // whole inside, else beside, else inside cut (range-labels.ts): an item labeled once
  lay.labeled.clear();
  labelsInside(lay, ctx.text, ctx.tokens, 0);
  if (lay.labelY >= 0) placeBeside(lay, ctx);
  labelsInside(lay, ctx.text, ctx.tokens, ctx.config.labels.minChars);
}

/** Labels beside, in the track's label lane: each item with no label inside, by priority; dropped when it cannot. */
function placeBeside(lay: RangesLay, ctx: LayoutContext) {
  const left = lay.left;
  const labeled = lay.labeled;
  const n = lay.items;
  const lane = ctx.labelLane ?? OwnLane;
  if (!ctx.labelLane) lane.reset();
  // the layer's marks, for the leaders of this and the later layers (a leader never crosses a mark above its item)
  for (let j = 0; j < lay.hn; j++) if (lay.hitem[j] >= 0) lane.addMark(lay.hx0[j], lay.hx1[j], lay.hlane[j]);
  // the candidates: an item once (its first row in view with a mark), by priority (the longer item first)
  lay.cand = grown(lay.cand, n);
  const cand = lay.cand;
  let c = 0;
  for (let i = 0; i < n; i++) {
    if (labeled.has(lay.keys[i]) || !lay.labels[i] || !(left[i] < Number.POSITIVE_INFINITY)) continue;
    labeled.add(lay.keys[i]);
    cand[c++] = i;
  }
  if (!c) return;
  // a run of one label (neighbours less than the label's width apart) is labeled once, at its start
  c = firstOfRuns(lay, ctx, c);
  SortStarts = lay.starts;
  SortEnds = lay.ends;
  sortArrayRange(cand, 0, c, longerFirst);
  SortStarts = SortEnds = NoNumbers;
  const x0 = lay.box.x;
  const x1 = lay.box.x + lay.box.width;
  for (let k = 0; k < c; k++) {
    const i = cand[k];
    const w = ctx.text.measure(lay.font, lay.labels[i] as string);
    let want = Math.max(left[i], x0);
    // the first free place from its start rightward, past the labels already placed in the track's lane
    let x = lane.free(want, w, LabelGapPx);
    // moved: its leader rises from its start, or (an item in a lower lane) from the first x of it with no mark above
    if (x - want >= 1 && lay.lanes[i] > 0)
      x = lane.free((want = lane.clearX(want, Math.min(lay.right[i], x1), lay.lanes[i])), w, LabelGapPx);
    // dropped past the box's end, or moved further than its own width (a leader is a short one), or when its leader
    // would cross a label or a mark (then it may go inside cut)
    const no = !(want < x1) || x + w > x1 || x - want > Math.max(w, 24);
    if (no || (x - want >= 1 && lane.crosses(want + 0.5, Math.max(0, lay.lanes[i]), x - 2))) {
      labeled.delete(lay.keys[i]);
      continue;
    }
    lane.add(x, w);
    const q = lay.ln++;
    if (q >= lay.lx.length) {
      lay.lx = grown(lay.lx, q + 1);
      lay.lw = grown(lay.lw, q + 1);
      lay.lwant = grown(lay.lwant, q + 1);
      lay.litem = grown(lay.litem, q + 1);
    }
    lay.lx[q] = x;
    lay.lw[q] = w;
    lay.lwant[q] = want;
    lay.litem[q] = i;
  }
}

const NoNumbers: readonly number[] = [];

/** The candidates' order by their first mark's x: a comparator of the module, reading this. */
let SortLeft: ArrayLike<number> = NoNumbers;
function leftFirst(o: ArrayLike<number>, a: number, b: number): number {
  return SortLeft[o[a]] - SortLeft[o[b]] || o[a] - o[b];
}

/**
 * The candidates `lay.cand[0..c)` without the followers in a run: an item whose label is the one before it in x, less
 * than a label's width from that item's end, is part of its run and goes unlabeled (its hover still names it). Returns
 * the count kept (the run's first items, not in order).
 */
function firstOfRuns(lay: RangesLay, ctx: LayoutContext, c: number): number {
  const cand = lay.cand;
  const left = lay.left;
  SortLeft = left;
  sortArrayRange(cand, 0, c, leftFirst);
  SortLeft = NoNumbers;
  const ends = lay.runEnds;
  ends.clear();
  let k = 0;
  for (let q = 0; q < c; q++) {
    const i = cand[q];
    const label = lay.labels[i] as string;
    const right = lay.right[i];
    const prev = ends.get(label);
    const w = ctx.text.measure(lay.font, label);
    if (prev !== undefined && left[i] - prev < w) {
      if (right > prev) ends.set(label, right);
      continue;
    }
    ends.set(label, right);
    cand[k++] = i;
  }
  return k;
}
