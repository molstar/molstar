/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The store's keys and tiles (_spec/README.md "Data", Store): a data entry is (source, system, level, tile), a node
// entry (source, key), each one string for `KeyedState`. Keys are written, never parsed: an entry keeps its parts.

/** Bins per tile, by default: a region at a level snaps to whole tiles, so overlapping requests share them. */
export const DefaultTileBins = 1024;

const Special = /[\\|#]/;

/** `|` and `#` separate the parts, so a part holding one (or `\`) escapes it: no two keys of other parts are equal. */
function part(s: string): string {
  return Special.test(s) ? s.replace(/[\\|#]/g, '\\$&') : s;
}

/** `<source>|<system>|<level>|<tile>`. */
export function dataKey(source: string, system: string, level: number, tile: number): string {
  return `${part(source)}|${part(system)}|${level}|${tile}`;
}

/** `<source>#<key>`. */
export function nodeKey(source: string, key: string): string {
  return `${part(source)}#${part(key)}`;
}

/** Positions one tile covers at `level`. */
export function tileWidth(level: number, tileBins = DefaultTileBins): number {
  return level * tileBins;
}

/** The tiles `[start, end)` touches at `level`: `[first, last]`, empty (`last < first`) for an empty range. */
export function tileSpan(start: number, end: number, level: number, tileBins = DefaultTileBins): [number, number] {
  const w = tileWidth(level, tileBins);
  const first = Math.floor(Math.max(0, start) / w);
  if (end <= start) return [first, first - 1];
  return [first, Math.ceil(end / w) - 1];
}

/** The positions of `tile` at `level`, clipped to the system's length when known. */
export function tileRegion(
  tile: number,
  level: number,
  tileBins = DefaultTileBins,
  length?: number,
): { start: number; end: number } {
  const w = tileWidth(level, tileBins);
  const start = tile * w;
  const end = length === undefined ? start + w : Math.min(start + w, length);
  return { start, end };
}

/**
 * The level that fits: the coarsest listed bin width no wider than the positions a device pixel covers; the finest
 * when every listed one is wider. `levels` ascending.
 */
export function fitLevel(levels: readonly number[], positionsPerPixel: number): number {
  let level = levels[0];
  for (const l of levels) {
    if (l <= positionsPerPixel) level = l;
    else break;
  }
  return level;
}
