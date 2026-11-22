/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { sortArrayRange } from '@molstar/core/data/util/sort';
import { projectRanges } from '@molstar/seqx-core/resolve';
import { boundedLevel, mergeNear, readTiles } from './layer-data.js';
import { collectItems, newItemRows } from './items.js';
import { noteLanes } from './known.js';
import type { LayoutContext, StoreReader } from './layer-type.js';
import type { LayerSource } from './layer-data.js';
import type { RangesPrep } from './ranges-lay.js';

// The lanes of a source's layer over its reachable region (_spec/README.md "Render", Track heights): the items of the
// section's whole extent, its tiles at the level that fits it, packed once they are all read; a track of a source is
// as tall as they need from the first settle, and a pan or a zoom keeps it.

/** Where a prep's reachable read stands: unread (absent), asked and waiting, or done. */
const Reach = new WeakMap<RangesPrep, { done: boolean }>();

/** The gap a gene model's label leaves before it, CSS px (`ranges-gene.ts`: its label and the gap to its model). */
const LabelGap = 4;

/**
 * Wants the tiles of the section's whole extent at the fit's level (every frame: a want set is the whole set) and, once
 * every one is read as items and the view has settled, notes the lanes they pack into (`noteLanes`). Returns whether the reachable lanes are known (the layer then notes none of its
 * own: its view's lanes are within them). A source whose fit level holds no items (a density's coarse tiles), one with
 * no known length or a failed tile is not known: the view's lanes are noted as before.
 */
export function noteReachable(prep: RangesPrep, src: LayerSource, ctx: LayoutContext, data: StoreReader): boolean {
  const kept = Reach.get(prep);
  if (kept?.done) return true;
  if (data.coverage(src.source) === undefined) return false;
  // the section's own pieces: its window of the system (a range's), never the whole system's length
  const whole = Array.from(ctx.axis.pieces);
  // the whole extent on the data's system, its stretches near one another read as one
  const regions = mergeNear(HalfOpenRanges.normalize(src.back ? projectRanges(src.back, whole) : whole));
  let span = 0;
  for (let i = 0; i + 1 < regions.length; i += 2) span += regions[i + 1] - regions[i];
  if (!(span > 0)) return false;
  const levels = data.coverage(src.source)?.levels ?? [1];
  const level = boundedLevel(levels, span / Math.max(1, ctx.box.width * ctx.axis.dpr), span, data);
  // the finest level holds the items: a coarser one is a density, whose lanes nobody knows
  if (level !== levels[0]) return false;
  const read = readTiles(ctx, data, src, regions, level, {
    layer: `${ctx.key}/reach`,
    priority: 'margin',
  });
  if (read.tooMany || read.failed.length) return false;
  // waiting for the rest: the layer's own lanes are not noted meanwhile (the track grows once, to the whole)
  if (read.pending.length || !ctx.settled) return true;
  const rows = newItemRows();
  const densities = collectItems(rows, read.entries, [], src.filter, src.base);
  if (densities.length) return false;
  noteLanes(prep, packed(prep, rows, ctx, span / Math.max(1, ctx.box.width)));
  Reach.set(prep, { done: true });
  return true;
}

/** The lanes first-fit by start packs the items into (a gene line's: an item's rows by its hull, a label's room before it). */
function packed(prep: RangesPrep, rows: ReturnType<typeof newItemRows>, ctx: LayoutContext, perPx: number): number {
  const gene = prep.line === 'gene';
  const room = gene && (ctx.labels === undefined || ctx.labels === 'inside');
  const starts: number[] = [];
  const ends: number[] = [];
  const labels: (string | undefined)[] = [];
  if (gene) {
    const at = new Map<string | number, number>();
    for (let j = 0; j < rows.n; j++) {
      const k = at.get(rows.keys[j]);
      if (k === undefined) {
        at.set(rows.keys[j], starts.length);
        starts.push(rows.starts[j]);
        ends.push(rows.ends[j]);
        labels.push(rows.labels[j]);
      } else {
        starts[k] = Math.min(starts[k], rows.starts[j]);
        ends[k] = Math.max(ends[k], rows.ends[j]);
      }
    }
  } else
    for (let j = 0; j < rows.n; j++) {
      starts.push(rows.starts[j]);
      ends.push(rows.ends[j]);
    }
  const n = starts.length;
  const font = ctx.config.itemFont;
  const from = new Float64Array(n);
  const to = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const label = room ? labels[i] : undefined;
    from[i] = starts[i] - (label ? (ctx.text.measure(font, label) + LabelGap) * perPx : 0);
    to[i] = Math.max(ends[i], starts[i] + 1);
  }
  const order = new Int32Array(n);
  for (let i = 0; i < n; i++) order[i] = i;
  sortArrayRange(order, 0, n, (o, a, b) => from[o[a]] - from[o[b]] || o[a] - o[b]);
  const lanes: number[] = [];
  for (let k = 0; k < n; k++) {
    const i = order[k];
    let lane = 0;
    while (lane < lanes.length && lanes[lane] > from[i]) lane++;
    lanes[lane] = to[i];
  }
  return Math.min(Math.max(1, lanes.length), ctx.config.ranges.maxLanes);
}
