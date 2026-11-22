/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A liftover (_spec/README.md "Render", A liftover; EC-212): data over one assembly drawn on a sequence of another
// through a map read by region (`map.source`). The map's runs are read on the section's sequence by the end there
// (`RegionRequest.side`: `b`, the map's `to`), the section's window carried back along them onto each far sequence,
// the data read there, and each row carried forward along the runs: a part that maps nowhere is not drawn (EC-21), a
// part with several images is drawn at each (EC-20). Loaded with the genome chunk: a liftover is over assemblies.

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import type { DataEntry, PairsList } from '@molstar/seqx-core/data';
import { type Run, RunsLinear } from '@molstar/seqx-core/resolve';
import {
  dataRegions,
  type LayerSource,
  type LiftSource,
  mergeNear,
  pairsLevels,
  placeholderSpans,
  readTiles,
} from './layer-data.js';
import { collectItems, type DensityTile, type ItemRows, type ItemsBase, newItemRows } from './items.js';
import type { LayoutContext, StoreReader } from './layer-type.js';

type Loaded = Extract<DataEntry, { status: 'loaded' }>;

/** A map's runs read on one system: each far system's (by the name the source gives it) as runs from it to the near. */
export interface LiftRuns {
  readonly byFar: Map<string, RunsLinear>;
  /** Near positions whose runs are not read yet, and those whose read failed. */
  readonly pending: number[];
  readonly failed: number[];
  /** Rows whose ends differ in length: no run (left out). */
  unequal: number;
}

/** The items level of a read of `span` positions: the finest whose tiles a want's bound holds. */
function itemsLevel(data: StoreReader, source: string, span: number): number {
  return pairsLevels(data.coverage(source)?.levels ?? [1], 1, span, data).items;
}

/**
 * The runs of `lift`'s map with their near end in `regions` (sorted `[start, end, …]` on the near system), read as
 * items tiles by that end: a row a run when its ends are of one length (else counted), `from` its far start, its
 * direction the ends' relative strand. Each row once across tiles (by its `row`).
 */
export function liftRuns(
  ctx: LayoutContext,
  data: StoreReader,
  lift: LiftSource,
  regions: readonly number[],
): LiftRuns {
  const near = lift.near;
  let span = 0;
  for (let i = 0; i + 1 < regions.length; i += 2) span += regions[i + 1] - regions[i];
  const read = readTiles(ctx, data, near, regions, itemsLevel(data, near.source, span), { items: true });
  const byName = new Map<string, Run[]>();
  const seen = new Set<number>();
  let unequal = 0;
  const onB = lift.side === 'b';
  for (const e of read.entries as Loaded[]) {
    const t = e.tables.pairs;
    if (t?.form !== 'list') continue;
    unequal += runsOf(t, near.system, onB, seen, byName);
  }
  const byFar = new Map<string, RunsLinear>();
  for (const [name, runs] of byName) byFar.set(name, new RunsLinear(runs));
  return { byFar, pending: read.pending, failed: read.failed, unequal };
}

/** A pairs list's rows as runs from their far end to their near end on `system` into `out`; the unequal rows' count. */
function runsOf(t: PairsList, system: string, onB: boolean, seen: Set<number>, out: Map<string, Run[]>): number {
  let unequal = 0;
  for (let i = 0; i < t.aStart.length; i++) {
    const row = t.row?.[i];
    if (row !== undefined) {
      if (seen.has(row)) continue;
      seen.add(row);
    }
    const aName = t.aSystem?.get(i) ?? system;
    const bName = t.bSystem?.get(i) ?? system;
    if ((onB ? bName : aName) !== system) continue;
    const length = t.aEnd[i] - t.aStart[i];
    if (length <= 0 || t.bEnd[i] - t.bStart[i] !== length) {
      unequal++;
      continue;
    }
    const far = onB ? aName : bName;
    const farStart = onB ? t.aStart[i] : t.bStart[i];
    const nearStart = onB ? t.bStart[i] : t.aStart[i];
    const dir: 1 | -1 = (t.aStrand?.[i] || 1) * (t.bStrand?.[i] || 1) < 0 ? -1 : 1;
    const run: Run = { from: farStart, to: dir === 1 ? nearStart : nearStart + length - 1, length, dir };
    let list = out.get(far ?? '');
    if (!list) out.set(far ?? '', (list = []));
    list.push(run);
  }
  return unequal;
}

/** What a lifted layer keeps between layouts: the far rows read, reused. */
const Scratch = new WeakMap<ItemsBase, { rows: ItemRows; buf: number[]; occ: Map<string | number, number> }>();

/**
 * A lifted layer's items into `lay.rows`, on the section's system (`readItems`' for `src.lift`): each far row carried
 * along the runs, a row an image (EC-20; a part with none is not drawn, EC-21), the placeholders where the runs or the
 * far data are not read. No density tables: the far data is read as items.
 */
export function readLifted(
  lay: ItemsBase,
  ctx: LayoutContext,
  data: StoreReader,
  src: LayerSource,
  fields: readonly string[],
): DensityTile[] {
  const lift = src.lift!;
  const axis = ctx.axis;
  const regions = dataRegions(src, axis, true);
  const runs = liftRuns(ctx, data, lift, regions);
  placeholderSpans(src, axis, runs.pending, lay.pending);
  placeholderSpans(src, axis, runs.failed, lay.failed);
  let s = Scratch.get(lay);
  if (!s) Scratch.set(lay, (s = { rows: newItemRows(), buf: [], occ: new Map() }));
  const rows = lay.rows;
  rows.n = 0;
  rows.filtered = 0;
  rows.clipped = 0;
  rows.leftOut = 0;
  rows.clamped = 0;
  s.occ.clear();
  const { buf, occ } = s;
  const back: number[] = [];
  for (const [name, linear] of runs.byFar) {
    const far = lift.far.get(name);
    if (!far) continue;
    // the near window carried back onto the far system
    back.length = 0;
    const inv = linear.inverse();
    for (let i = 0; i + 1 < regions.length; i += 2) inv.targets(regions[i], regions[i + 1], back);
    const farRegions = mergeNear(HalfOpenRanges.normalize(back));
    let span = 0;
    for (let i = 0; i + 1 < farRegions.length; i += 2) span += farRegions[i + 1] - farRegions[i];
    const read = readTiles(ctx, data, far, farRegions, itemsLevel(data, far.source, span), { items: true });
    for (const list of [read.pending, read.failed]) {
      buf.length = 0;
      for (let i = 0; i + 1 < list.length; i += 2) linear.targets(list[i], list[i + 1], buf);
      placeholderSpans(src, axis, buf, list === read.pending ? lay.pending : lay.failed);
    }
    const wrap = far.length !== undefined && far.wraps ? { length: far.length, cyclic: false } : undefined;
    collectItems(s.rows, read.entries, fields, far.filter, far.base, wrap);
    carry(s.rows, linear, rows, buf, occ);
    rows.filtered += s.rows.filtered;
    rows.clipped += s.rows.clipped;
    rows.leftOut += s.rows.leftOut;
  }
  rows.link.length = rows.n;
  rows.link.fill(-1);
  return [];
}

/** Each far row's images along `linear` appended to `out` as rows of its item (an image a row, in source order). */
function carry(from: ItemRows, linear: RunsLinear, out: ItemRows, buf: number[], occ: Map<string | number, number>) {
  for (let j = 0; j < from.n; j++) {
    buf.length = 0;
    linear.targets(from.starts[j], from.ends[j], buf);
    if (!buf.length) continue;
    const pieces = HalfOpenRanges.normalize(buf);
    const key = from.keys[j];
    for (let k = 0; k + 1 < pieces.length; k += 2) {
      const n = out.n++;
      const o = occ.get(key) ?? 0;
      occ.set(key, o + 1);
      out.keys[n] = key;
      out.starts[n] = pieces[k];
      out.ends[n] = pieces[k + 1];
      out.values[n] = from.values[j];
      out.labels[n] = from.labels[j];
      out.occ[n] = o;
      out.fields[n] = from.fields[j];
      out.row[n] = from.row[j];
      out.itemStarts[n] = -1;
      out.itemEnds[n] = -1;
    }
  }
}
