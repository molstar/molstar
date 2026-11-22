/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesDensity, RangesList, RegionRequest } from '../../data/protocol.js';

// The levels of a tab-line source read by region (_spec/README.md "Data", Genome runtimes' internals: levels): listed
// per source from its index's bytes (a level only when every tile of it on every system fits `maxRegionBytes`, the cost
// read from the index, never the data), and above level 1 a density counted exactly from the tile's own records. Shared
// by the BED, GFF and VCF runtimes, not published (`./runtimes/internal/*`).

/** The candidate levels: tile widths of 1, 16, 256 and 4,096 kb at 1,024 bins a tile. */
export const LevelCandidates: readonly number[] = [1, 16, 256, 4096];

/** What bounds a source's levels: the coarsest listed, the store's bins a tile, a tile's compressed bytes. */
export interface LevelBounds {
  /** The coarsest level listed (4,096: 4 Mb tiles). */
  readonly maxLevel: number;
  /** Bins a tile (the store's `tileBins`, 1,024): a level's tile is `level × tileBins` positions. */
  readonly tileBins: number;
  /** A tile's chunks, compressed, together (io's `maxRegionBytes`, 8 MiB). */
  readonly maxRegionBytes: number;
}

/** A tile's cost from the index (`regionBytes`), nothing read. */
export type TileCost = (name: string, start: number, end: number) => number;

/** The longest span an index covers (a `.tbi`'s 2^29), never past int32. */
export function indexSpan(minShift: number, depth: number): number {
  return Math.min(2 ** (minShift + 3 * depth), 0x7fffffff);
}

/**
 * Where a system's indexed lines end, at most: the end of the last tile of `width` with a cost (its length unknown: no
 * `##contig`, no `##sequence-region`), scanning the index's span at that width (a `.tbi`'s 2^29 at 4 Mb: 128 tiles).
 */
export function indexedExtent(cost: TileCost, name: string, width: number, span: number): number {
  let extent = 0;
  for (let s = 0; s < span; s += width) if (cost(name, s, Math.min(span, s + width)) > 0) extent = s + width;
  return Math.min(extent, span);
}

/**
 * The levels a source lists: 1, and each candidate up to `maxLevel` whose every tile (`level × tileBins` positions from
 * 0, on every system to its length) costs at most `maxRegionBytes`. Each level is checked tile by tile, the first tile
 * past the bound ending its check; level 1 is always listed (a tile past the bound is answered `too-large`). Why every
 * tile and from the index: the bytes a tile costs are known exactly before a read, its rows are not (EC-171).
 */
export function tabixLevels(
  cost: TileCost,
  systems: readonly { readonly name: string; readonly extent: number }[],
  bounds: LevelBounds,
): number[] {
  const out = [1];
  for (const level of LevelCandidates) {
    if (level === 1 || level > bounds.maxLevel) continue;
    const width = level * bounds.tileBins;
    let fits = true;
    for (const { name, extent } of systems) {
      for (let s = 0; fits && s < extent; s += width)
        fits = cost(name, s, Math.min(extent, s + width)) <= bounds.maxRegionBytes;
      if (!fits) break;
    }
    if (fits) out.push(level);
  }
  return out;
}

/**
 * A ranges list's density at `level` over `[start, end)` (a tile): each item (by `id`, else `row`, else its row) counted
 * once in every bin its extent touches (`itemStart`, `itemEnd` when the rows are pieces, else the hull of its rows), the
 * bins inside the region only; bins with no item left out. Counted exactly from the tile's records, never estimated.
 * One pass over the rows and one over the bins (counts by difference); with `unique` (a row an item: VCF's alleles) no
 * identity is looked up, so nothing is allocated a record.
 */
export function densityOf(
  list: RangesList,
  request: Pick<RegionRequest, 'start' | 'end' | 'level'>,
  unique = false,
): RangesDensity {
  const { level } = request;
  const lo = Math.floor(request.start / level);
  const hi = Math.max(lo, Math.ceil(request.end / level));
  const diff = new Float64Array(hi - lo + 1);
  const n = list.start.length;
  const add = (s: number, e: number) => {
    const a = Math.max(lo, Math.floor(s / level));
    const b = Math.min(hi - 1, Math.floor((Math.max(e, s + 1) - 1) / level));
    if (b < a) return;
    diff[a - lo]++;
    diff[b - lo + 1]--;
  };
  const ids = list.id;
  const rows = list.row;
  const is = list.itemStart;
  const ie = list.itemEnd;
  // an identity seen twice (pieces of one item): its hull, counted once after the pass
  let hulls: Map<string | number, [number, number]> | undefined;
  let first: Map<string | number, number> | undefined;
  for (let i = 0; i < n; i++) {
    const s = is ? is[i] : list.start[i];
    const e = ie ? ie[i] : list.end[i];
    if (unique) {
      add(s, e);
      continue;
    }
    const id = ids?.get(i);
    const key: string | number | undefined = id !== null && id !== undefined ? id : rows ? rows[i] : undefined;
    if (key === undefined) {
      add(s, e);
      continue;
    }
    first ??= new Map();
    const at = first.get(key);
    if (at === undefined) {
      first.set(key, i);
      continue;
    }
    hulls ??= new Map();
    const h = hulls.get(key);
    if (h) {
      h[0] = Math.min(h[0], s);
      h[1] = Math.max(h[1], e);
    } else {
      const s0 = is ? is[at] : list.start[at];
      const e0 = ie ? ie[at] : list.end[at];
      hulls.set(key, [Math.min(s0, s), Math.max(e0, e)]);
    }
  }
  if (first)
    for (const [key, i] of first) {
      const h = hulls?.get(key);
      if (h) add(h[0], h[1]);
      else add(is ? is[i] : list.start[i], ie ? ie[i] : list.end[i]);
    }
  let bins = 0;
  for (let k = 0, run = 0; k < hi - lo; k++) {
    run += diff[k];
    if (run > 0) bins++;
  }
  const bin = new Int32Array(bins);
  const count = new Float32Array(bins);
  for (let k = 0, run = 0, j = 0; k < hi - lo; k++) {
    run += diff[k];
    if (run > 0) {
      bin[j] = lo + k;
      count[j++] = run;
    }
  }
  return { form: 'density', bin, count };
}
