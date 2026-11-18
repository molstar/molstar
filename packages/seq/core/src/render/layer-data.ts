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
}

/** The store's source of a placed layer's data; undefined when its node has no access plan. */
export function layerSource(model: StaticModel, layer: ResolvedLayer): LayerSource | undefined {
  const data = layer.data;
  if (!data || !layer.placed) return undefined;
  // a template instance's data is the overlay's: its plan and length travel with it
  let access = data.access ?? model.systems.get(data.address)?.access;
  if (!access) {
    let a = model.annotations.get(data.address);
    // a query reads the annotation it starts at (its filter is a later slice's)
    if (!a) for (const [address, entry] of model.annotations) if (data.address.startsWith(`${address}/`)) a = entry;
    access = a?.access;
  }
  if (!access) return undefined;
  const plan =
    access.kind === 'inline'
      ? { source: access.source.address, system: access.system, shape: access.shape }
      : { source: access.source, system: access.system, shape: access.shape };
  const path = layer.path;
  const length = data.length ?? model.systems.get(data.system)?.length;
  return {
    ...plan,
    dataSystem: data.system,
    ...(length !== undefined ? { length } : {}),
    ...(path ? { path, back: inversePath(path) } : {}),
    ratio: path?.step === 1 ? 1 / 3 : path?.step === -1 ? 3 : 1,
  };
}

/** The data's positions in view (`[start, end, …]`, sorted, merged): the section's carried back along the path. */
export function dataRegions(src: LayerSource, axis: AxisTransform): number[] {
  const vis = axis.visible([]);
  if (!src.back) return Array.from(HalfOpenRanges.normalize(vis));
  return Array.from(HalfOpenRanges.normalize(projectRanges(src.back, vis)));
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

/** The level that fits the view on the data's system: the coarsest bin no wider than a device pixel's positions. */
export function layerLevel(src: LayerSource, data: StoreReader, axis: AxisTransform, finest = false): number {
  const levels = data.coverage(src.source)?.levels ?? [1];
  return finest ? levels[0] : fitLevel(levels, src.ratio / axis.devicePx);
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
