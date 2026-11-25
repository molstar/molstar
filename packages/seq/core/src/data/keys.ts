/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The store's keys and tiles (_spec/README.md "Data", Store): a data entry is (source, system, level, tile), a 2D one
// with its second axis (system, tile), a node entry (source, key), each one string for `KeyedState`. Keys are written,
// never parsed: an entry keeps its parts.

/** Bins per tile, by default: a region at a level snaps to whole tiles, so overlapping requests share them. */
export const DefaultTileBins = 1024;

/** Bins a side of a 2D (grid) tile, by default: a dense tile of 256² float cells is 256 kB, one image tile. */
export const DefaultGridTileBins = 256;

const Special = /[\\|#]/;
const SystemSpecial = /[\\|#:]/;

/** `|` and `#` separate the parts, so a part holding one (or `\`) escapes it: no two keys of other parts are equal. */
function part(s: string): string {
  return Special.test(s) ? s.replace(/[\\|#]/g, '\\$&') : s;
}

/** A system part: `:` escaped too, so a sided one (`a:chr1`) never equals an unsided name. */
function systemPart(s: string, side: 'a' | 'b' | undefined): string {
  const p = SystemSpecial.test(s) ? s.replace(/[\\|#:]/g, '\\$&') : s;
  return side ? `${side}:${p}` : p;
}

/** What a data key carries besides its source, system, level and tile (_spec/README.md "Data", Store, Keys). */
export interface KeyParts {
  /** A sided pairs source's side: the system part `a:chr1`. */
  readonly side?: 'a' | 'b';
  /** An items tile: the level part `i<level>`, never a bins tile's key at the same level. */
  readonly items?: boolean;
  /** A 2D tile's second axis: two more parts, so a 2D key never equals a 1D one. */
  readonly y?: { readonly system: string; readonly tile: number };
}

/** `<source>|<system>|<level>|<tile>`, and `|<y system>|<y tile>` for a 2D tile. */
export function dataKey(source: string, system: string, level: number, tile: number, parts?: KeyParts): string {
  const key = `${part(source)}|${systemPart(system, parts?.side)}|${parts?.items ? 'i' : ''}${level}|${tile}`;
  return parts?.y ? `${key}|${systemPart(parts.y.system, undefined)}|${parts.y.tile}` : key;
}

/**
 * A symmetric grid's tile folded to the one served (_spec/README.md "Data", Store): the tile at or above the diagonal,
 * its x tile not after its y tile by `(system order, tile)`; `transposed` when the asked tile is its mirror.
 */
export function foldTile<T extends { readonly system: string; readonly tile: number }>(
  x: T,
  y: T,
  order: (system: string) => number,
): { x: T; y: T; transposed: boolean } {
  const ox = order(x.system);
  const oy = order(y.system);
  const after = ox > oy || (ox === oy && x.tile > y.tile);
  return after ? { x: y, y: x, transposed: true } : { x, y, transposed: false };
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
