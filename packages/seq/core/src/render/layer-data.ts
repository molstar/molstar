/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { type DataEntry, fitLevel, type Shape } from '@molstar/seqx-core/data';
import { inversePath, type Path, projectRanges, type StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { compileRowFilter, type RowFilter } from './filter.js';
import type { LayoutContext, StoreReader } from './layer-type.js';

// How a layer reads its data (_spec/README.md "Render", Data): the access plan of its data's node names the store's
// source and system; the region in view is carried back along the layer's path onto the data's system; the level fits
// the positions a device pixel covers there; the tiles' entries are read, and those not loaded are placeholders.

/** Where a layer's data comes from in the store. */
export interface LayerSource {
  /** The store's source (an inline node's own address, or a source node's). */
  readonly source: string;
  /** The name the source knows the system by. */
  readonly system: string;
  /** The data's system in the spec. */
  readonly dataSystem: string;
  readonly shape: Shape;
  /** The data's system's length, when known. */
  readonly length?: number;
  /** The path from the data's system to the section's (none: the same) and back. */
  readonly path?: Path;
  readonly back?: Path;
  /** Data positions per section position along the path: 1/3 from residues onto nucleotides, 3 back, else 1. */
  readonly ratio: number;
  /** The filter the layer's query puts on its rows (an item layer's), compiled once; `unread`: a form it cannot read. */
  readonly filter?: RowFilter;
  readonly unreadFilter?: boolean;
}

/** The store's source of a placed layer's data; undefined when its node has no access plan. */
export function layerSource(model: StaticModel, layer: ResolvedLayer): LayerSource | undefined {
  const data = layer.data;
  if (!data || !layer.placed) return undefined;
  // a template instance's data is the overlay's: its plan and length travel with it
  let access = data.access ?? model.systems.get(data.address)?.access;
  if (!access) {
    let a = model.annotations.get(data.address);
    // a query reads the annotation it starts at; its form after it filters the rows (`filter`)
    if (!a) for (const [address, entry] of model.annotations) if (data.address.startsWith(`${address}/`)) a = entry;
    access = a?.access;
  }
  if (!access) return undefined;
  const plan =
    access.kind === 'inline'
      ? { source: access.source.address, system: access.system, shape: access.shape }
      : { source: access.source, system: access.system, shape: access.shape };
  const path = layer.path;
  const filter = data.filter ? compileRowFilter(data.filter) : undefined;
  const length = data.length ?? model.systems.get(data.system)?.length;
  return {
    ...plan,
    dataSystem: data.system,
    ...(length !== undefined ? { length } : {}),
    ...(path ? { path, back: inversePath(path) } : {}),
    ratio: path?.step === 1 ? 1 / 3 : path?.step === -1 ? 3 : 1,
    ...(filter ? { filter } : {}),
    ...(data.filter && !filter ? { unreadFilter: true } : {}),
  };
}

/**
 * The data's positions in view (`[start, end, …]`, sorted, merged): the section's carried back along the path; with
 * `hidden`, the stretches the breaks in view hide too (what a layer draws as their ticks). Near ones merged (`near`,
 * `mergeNear`: what a layer reads); without, exactly those in view (what a layer walks position by position).
 */
export function dataRegions(src: LayerSource, axis: AxisTransform, hidden = false, near = true): number[] {
  const vis = axis.visible([]);
  if (hidden) {
    const b = axis.breaks([]);
    for (let i = 0; i + 2 < b.length; i += 3)
      if (b[i + 2] >= b[i + 1]) vis.push(b[i + 1], b[i + 2]);
      else vis.push(b[i + 1], axis.length ?? b[i + 1], 0, b[i + 2]);
  }
  const exact = HalfOpenRanges.normalize(src.back ? projectRanges(src.back, vis) : vis);
  return near ? mergeNear(exact) : Array.from(exact);
}

/** Gaps of at most this many positions between regions are read with them (small: an intron of an exon pair stays apart). */
export const ReadGap = 16;

/**
 * Sorted, merged regions with the gaps of at most `ReadGap` positions between them closed: a collapsed axis's pieces
 * in view (150,000 single positions, EC-33) are read as a few regions, not a want each.
 */
export function mergeNear(ranges: ArrayLike<number>, gap = ReadGap): number[] {
  const out: number[] = [];
  for (let i = 0; i + 1 < ranges.length; i += 2) {
    const n = out.length;
    if (n && ranges[i] - out[n - 1] <= gap) out[n - 1] = Math.max(out[n - 1], ranges[i + 1]);
    else out.push(ranges[i], ranges[i + 1]);
  }
  return out;
}

/** A region on the data's system carried onto the section's as pieces (`[start, end, …]`), into `out`. */
export function toSection(src: LayerSource, s: number, e: number, out: number[]): number[] {
  out.length = 0;
  if (!src.path) {
    out.push(s, e);
    return out;
  }
  return projectRanges(src.path, [s, e], out);
}

/**
 * The level that fits the view on the data's system: the coarsest bin no wider than `coarser` device pixels' positions
 * (a values layer's `values.coarser`, else 1), and coarser still while the view's span at it would want more tiles
 * than the store's bound (`maxTilesPerWant`, where the reader says it: a whole chromosome at a high device pixel ratio
 * reads its coarse level, never a refused one); with `finest`, the finest level.
 */
export function layerLevel(
  src: LayerSource,
  data: StoreReader,
  axis: AxisTransform,
  finest = false,
  coarser = 1,
): number {
  const levels = data.coverage(src.source)?.levels ?? [1];
  if (finest) return levels[0];
  let level = fitLevel(levels, (coarser * src.ratio) / axis.devicePx);
  const s = data.settings;
  if (s) {
    const span = (axis.u1 - axis.u0) * src.ratio;
    for (let i = levels.indexOf(level); i + 1 < levels.length && span / (s.tileBins * level) + 1 > s.maxTilesPerWant; )
      level = levels[++i];
  }
  return level;
}

/** The tiles' state over the regions read: loaded entries, and the data positions still pending or failed. */
export interface TileRead {
  readonly entries: Extract<DataEntry, { status: 'loaded' }>[];
  /** `[start, end, …]` on the data's system. */
  readonly pending: number[];
  readonly failed: number[];
}

/**
 * Wants the regions (the view at `view` priority, a margin each side) and reads the tiles of the view's. A tile not
 * asked for, pending or evicted is pending (EC-169: never drawn as empty); a failed one is failed; an unknown source
 * (not attached yet) is pending over the regions.
 */
export function readTiles(
  ctx: LayoutContext,
  data: StoreReader,
  src: LayerSource,
  regions: readonly number[],
  level: number,
): TileRead {
  const out: TileRead = { entries: [], pending: [], failed: [] };
  const known = data.coverage(src.source) !== undefined;
  // two regions can share a tile (a cyclic view wrapping a short system): each entry is read once
  const seen = new Set<string>();
  for (let i = 0; i + 1 < regions.length; i += 2) {
    const start = regions[i];
    const end = regions[i + 1];
    const want = { source: src.source, system: src.system, start, end, level, priority: 'view' as const };
    ctx.want(want);
    const m = Math.ceil((end - start) * ctx.config.margin);
    if (m > 0) {
      const max = src.length ?? Number.MAX_SAFE_INTEGER;
      if (start > 0) ctx.want({ ...want, start: Math.max(0, start - m), end: start, priority: 'margin' });
      if (end < max) ctx.want({ ...want, start: end, end: Math.min(max, end + m), priority: 'margin' });
    }
    if (!known) {
      out.pending.push(start, end);
      continue;
    }
    for (const key of data.keysOf(want)) {
      const e = data.get(key);
      if (e?.status === 'loaded') {
        if (!seen.has(key)) out.entries.push(e);
        seen.add(key);
      } else if (e?.status === 'empty') continue;
      else if (e?.status === 'failed') out.failed.push(Math.max(start, e.region.start), Math.min(end, e.region.end));
      else if (e) out.pending.push(Math.max(start, e.region.start), Math.min(end, e.region.end));
      else out.pending.push(start, end);
    }
  }
  return out;
}

/** Placeholders' x extents: data positions carried onto the section and to x (`x0, x1, cut` triples), into `out`. */
export function placeholderSpans(src: LayerSource, axis: AxisTransform, ranges: readonly number[], out: number[]) {
  const pieces: number[] = [];
  const merged = HalfOpenRanges.normalize(ranges);
  for (let i = 0; i + 1 < merged.length; i += 2) {
    toSection(src, merged[i], merged[i + 1], pieces);
    for (let k = 0; k + 1 < pieces.length; k += 2) axis.spans(pieces[k], pieces[k + 1], out);
  }
  return out;
}
