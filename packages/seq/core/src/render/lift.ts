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
import type { DataEntry, PairsList, Want } from '@molstar/seqx-core/data';
import {
  type AccessPlan,
  type AnnotationEntry,
  type Run,
  RunsLinear,
  type StaticModel,
} from '@molstar/seqx-core/resolve';
import type { ResolvedLayer, ResolvedView } from '@molstar/seqx-core/view';
import { regionMaps } from '@molstar/seqx-core/view/lift';
import { spanText } from './accessibility.js';
import { compileRowFilter } from './filter.js';
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

/** What a projection through a region map reads of the store: the tiles held now, never a want. */
export interface RegionReader {
  coverage(source: string): { readonly levels?: readonly number[] } | undefined;
  keysOf(want: Want): string[];
  get(key: string): DataEntry | undefined;
}

/** The side a sequence is on in a region map's annotation: `a` its `from` assembly's, `b` its `to`'s. */
function sideOf(model: StaticModel, m: AnnotationEntry, system: string): 'a' | 'b' | undefined {
  if (model.assemblies.get(m.target)?.systems.includes(system)) return 'a';
  if (m.ySystem !== undefined && model.assemblies.get(m.ySystem)?.systems.includes(system)) return 'b';
  return undefined;
}

/**
 * `ranges` on `from` carried onto `to` through a map read by region joining them (`RegionProjector.project`): the runs
 * with an end on `from` in the items tiles the store holds at any level (read by that end), each range's images through
 * the runs to `to` appended to `out` (EC-20: every image); false when no such map joins them.
 */
export function projectLifted(
  model: StaticModel,
  data: RegionReader,
  from: string,
  ranges: ArrayLike<number>,
  to: string,
  out: number[],
): boolean {
  for (const address of regionMaps(model.spec)) {
    const m = model.annotations.get(address);
    const near = m?.bySystem?.get(from);
    const far = m?.bySystem?.get(to);
    if (!m || !near || !far || near.kind !== 'source') continue;
    const side = sideOf(model, m, from);
    if (!side || side === sideOf(model, m, to)) continue;
    const byName = new Map<string, Run[]>();
    const seen = new Set<number>();
    for (const level of data.coverage(near.source)?.levels ?? [1])
      for (let i = 0; i + 1 < ranges.length; i += 2) {
        const want: Want = {
          source: near.source,
          system: near.system,
          start: ranges[i],
          end: ranges[i + 1],
          level,
          items: true,
          ...(near.side ? { side: near.side } : {}),
        };
        for (const key of data.keysOf(want)) {
          const e = data.get(key);
          const t = e?.status === 'loaded' ? e.tables.pairs : undefined;
          if (t?.form === 'list') runsOf(t, near.system, side === 'b', seen, byName);
        }
      }
    const runs = byName.get(far.system);
    if (!runs) return true;
    // the runs go from the far end to the near: read back, from `from` to `to`
    const back = new RunsLinear(runs).inverse();
    for (let i = 0; i + 1 < ranges.length; i += 2) back.targets(ranges[i], ranges[i + 1], out);
    return true;
  }
  return false;
}
/** The annotation a layer's data reads: its address's, else the one a query starts at. */
function annotationOf(model: StaticModel, address: string) {
  const a = model.annotations.get(address);
  if (a) return a;
  for (const [at, entry] of model.annotations) if (address.startsWith(`${at}/`)) return entry;
  return undefined;
}

/** A plan's store source and system (an inline node's own address). */
const planOf = (p: AccessPlan) =>
  p.kind === 'inline'
    ? { source: p.source.address, system: p.system, shape: p.shape }
    : { source: p.source, system: p.system, shape: p.shape, ...(p.side ? { side: p.side } : {}) };

/**
 * A lifted layer's source (_spec/README.md "Render", A liftover): its rows on the section's system, read as the map's
 * runs there and the data along them on each far system the map's source names and the data's annotation reads.
 */
export function liftSource(model: StaticModel, data: NonNullable<ResolvedLayer['data']>): LayerSource | undefined {
  const lift = data.lift!;
  const m = model.annotations.get(lift.map);
  const a = annotationOf(model, data.address);
  if (!m?.bySystem || !a?.bySystem) return undefined;
  const filter = data.filter ? compileRowFilter(data.filter) : undefined;
  const base = data.base ? compileRowFilter(data.base) : undefined;
  const length = model.systems.get(data.system)?.length;
  const len = (n: number | undefined) => (n !== undefined ? { length: n } : {});
  const far = new Map<string, LayerSource>();
  for (const [system, plan] of a.bySystem) {
    const named = m.bySystem.get(system);
    if (!named) continue;
    const n = model.systems.get(system)?.length;
    far.set(named.system, {
      ...planOf(plan),
      dataSystem: system,
      ...len(n),
      ...(n !== undefined && plan.kind !== 'inline' ? { wraps: true } : {}),
      ratio: 1,
      ...(filter ? { filter } : {}),
      ...(base ? { base } : {}),
    });
  }
  const near: LayerSource = { ...planOf(lift.near), dataSystem: data.system, ...len(length), ratio: 1 };
  return {
    source: near.source,
    system: near.system,
    dataSystem: data.system,
    shape: 'ranges',
    ...len(length),
    ratio: 1,
    ...(filter ? { filter } : {}),
    ...(base ? { base } : {}),
    ...((data.filter && !filter) || (data.base && !base) || data.unread ? { unreadFilter: true } : {}),
    lift: { map: lift.map, side: lift.side, near, far },
  };
}

/** The most places a card names: past it, the count alone. */
const MaxNamed = 4;

/**
 * `[start, end)` on `system` carried through the maps read by region onto every block `view` shows (`RegionProjector`'s
 * `places`): its places as `2 places: chr1 62,001–68,000, chr2 62,001–68,000` when there are at least two, else
 * undefined (one is the projected hover itself). A place is an image's pieces nearer than the item's length joined (a
 * chain's small gap splits no place).
 */
export function liftPlaces(
  model: StaticModel,
  data: RegionReader,
  view: ResolvedView,
  system: string,
  start: number,
  end: number,
): string | undefined {
  if (!regionMaps(model.spec).length) return undefined;
  const named: string[] = [];
  const seen = new Set<string>([system]);
  const buf: number[] = [];
  const range = [start, end];
  let n = 0;
  for (const s of view.sections)
    for (const to of s.blocks ? s.blocks.map((b) => b.axis.system) : [s.axis.system]) {
      if (seen.has(to)) continue;
      seen.add(to);
      buf.length = 0;
      if (!projectLifted(model, data, system, range, to, buf) || !buf.length) continue;
      const r = HalfOpenRanges.normalize(buf);
      const name = model.systems.get(to)?.name ?? to;
      for (let k = 0; k + 1 < r.length; n++) {
        const s0 = r[k];
        let e = r[k + 1];
        for (k += 2; k + 1 < r.length && r[k] - e <= end - start; k += 2) e = r[k + 1];
        if (named.length < MaxNamed) named.push(`${name} ${spanText(s0, e)}`);
      }
    }
  if (n < 2) return undefined;
  return `${n} places: ${named.join(', ')}${n > named.length ? ', …' : ''}`;
}
