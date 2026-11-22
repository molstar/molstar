/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { dataKey, fitLevel, foldTile, tileRegion, tileSpan } from './keys.js';
import { ItemShapes, type OpenedSource, type RegionRequest, type Shape } from './protocol.js';
import type { StoreSources } from './store-sources.js';
import { Rank, type StoreSettings, type TileRegion, type Want, type Wanter } from './store-types.js';
import { sourceNameOf } from './system-names.js';

// The data store's want sets and their planning (_spec/README.md "Data", Store): a consumer's wants made into the
// tiles they cover, each at its best priority. A 2D want's tiles are its x tiles × its y tiles, a symmetric grid's
// folded to the tile served; a want past `maxTilesPerWant` is not planned, nor a layer's wants past `maxTilesPerLayer`
// tiles or `maxBytesPerLayer` of read cost together (said by the caller).

/** A tile a want set needs, at its best rank. */
export interface PlannedTile {
  rank: number;
  readonly region: TileRegion;
  readonly shapes: readonly Shape[];
}

/** One tile of a want: its key, its region (the served one), and whether the asked tile is its mirror. */
export interface WantTile {
  readonly key: string;
  readonly region: TileRegion;
  readonly transposed: boolean;
}

/** A sided system's name in a source's names (`StoreSources`): its side, a separator no name holds, the name. */
export function sidedName(side: 'a' | 'b' | undefined, name: string): string {
  return side ? `${side}\u0001${name}` : name;
}

/** The tiles of a want: its systems' coverage names, the level read, the tiles' spans; undefined for an unknown one. */
export function wantTiles(w: Want, sources: StoreSources, settings: StoreSettings) {
  const s = sources.get(w.source);
  const system = s && sourceNameOf(s.names, sidedName(w.side, w.system));
  if (!s || system === undefined) return undefined;
  const cov = s.opened.coverage;
  const level = cov.levels.includes(w.level) ? w.level : fitLevel(cov.levels, w.level);
  const length = cov.systems.find((x) => x.name === system && x.side === w.side)?.length;
  const bins = w.y ? settings.gridTileBins : settings.tileBins;
  const end = length === undefined ? w.end : Math.min(w.end, length);
  const [first, last] = tileSpan(w.start, end, level, bins);
  const own = cov.shapes.filter((x) => x.system === system);
  if (!w.y) {
    // a 1D tile asks every shape of its system but grids (2D), an items tile its list shapes alone
    const shapes = own.map((x) => x.shape).filter((x) => x !== 'grid' && (!w.items || ItemShapes.includes(x)));
    return { system, level, length, first, last, shapes, bins, count: Math.max(0, last - first + 1) };
  }
  const ySystem = sourceNameOf(s.names, w.y.system);
  if (ySystem === undefined) return undefined;
  const grids = own.filter((x) => x.shape === 'grid' && (x.y === '*' || (x.y ?? [system]).includes(ySystem)));
  if (!grids.length) return undefined;
  const yLength = cov.systems.find((x) => x.name === ySystem)?.length;
  const yEnd = yLength === undefined ? w.y.end : Math.min(w.y.end, yLength);
  const [yFirst, yLast] = tileSpan(w.y.start, yEnd, level, bins);
  const order = (name: string) => cov.systems.findIndex((x) => x.name === name);
  const y = { system: ySystem, length: yLength, first: yFirst, last: yLast, symmetric: grids.some((g) => g.symmetric) };
  const count = Math.max(0, last - first + 1) * Math.max(0, yLast - yFirst + 1);
  return { system, level, length, first, last, shapes: ['grid'] as Shape[], bins, count, y, order };
}

type Tiles = NonNullable<ReturnType<typeof wantTiles>>;

/** The tiles of a planned want, x-major; a symmetric grid's tiles folded (each served tile once, its first mirror). */
export function eachTile(w: Want, t: Tiles): WantTile[] {
  const out: WantTile[] = [];
  const side = w.side;
  const items = w.items === true ? true : undefined;
  for (let tile = t.first; tile <= t.last; tile++) {
    const { start, end } = tileRegion(tile, t.level, t.bins, t.length);
    if (!t.y) {
      const key = dataKey(w.source, t.system, t.level, tile, { side, items });
      const region = { source: w.source, system: t.system, level: t.level, tile, start, end, side, items };
      out.push({ key, region: strip(region), transposed: false });
      continue;
    }
    for (let yt = t.y.first; yt <= t.y.last; yt++) {
      const x = { system: t.system, tile, length: t.length };
      const y = { system: t.y.system, tile: yt, length: t.y.length };
      const f = t.y.symmetric ? foldTile(x, y, t.order!) : { x, y, transposed: false };
      const xr = tileRegion(f.x.tile, t.level, t.bins, f.x.length);
      const yr = tileRegion(f.y.tile, t.level, t.bins, f.y.length);
      const key = dataKey(w.source, f.x.system, t.level, f.x.tile, { y: f.y });
      const region: TileRegion = {
        source: w.source,
        system: f.x.system,
        level: t.level,
        tile: f.x.tile,
        start: xr.start,
        end: xr.end,
        y: { system: f.y.system, tile: f.y.tile, start: yr.start, end: yr.end },
      };
      out.push({ key, region, transposed: f.transposed });
    }
  }
  return out;
}

/** A region without its absent optional parts (an entry's region compares by its own fields only). */
function strip(r: TileRegion): TileRegion {
  const { side, items, ...rest } = r;
  return { ...rest, ...(side && { side }), ...(items && { items }) };
}

/** The request asking a tile for its shapes. */
export function tileRequest(region: TileRegion, shapes: readonly Shape[]): RegionRequest {
  const { system, start, end, level, y, side, items } = region;
  return {
    system,
    start,
    end,
    level,
    shapes,
    ...(y && { y: { system: y.system, start: y.start, end: y.end } }),
    ...(side && { side }),
    ...(items && { items }),
  };
}

/** Why a want, or its layer's wants together, are not planned: past which bound, with its numbers. */
export interface TooMany {
  readonly want: Want;
  readonly bound: 'want' | 'layer' | 'cost';
  readonly value: number;
  readonly max: number;
}

/** What a layer's wants cover together: tiles, and the read cost of those whose source knows it. */
export interface LayerLoad {
  tiles: number;
  cost: number;
}

/** The read cost of a tile from its source (0 when it knows none). */
function costOf(opened: OpenedSource, region: TileRegion, shapes: readonly Shape[]): number {
  const c = opened.readCost?.(tileRequest(region, shapes));
  return typeof c === 'number' && c > 0 ? c : 0;
}

/**
 * A consumer's wants as tiles by key; `tooMany` is told of each want past `maxTilesPerWant`, and of each want of a
 * layer whose wants together pass `maxTilesPerLayer` tiles or `maxBytesPerLayer` of read cost (none of them planned).
 */
export function planWants(
  wants: readonly Want[],
  sources: StoreSources,
  settings: StoreSettings,
  tooMany: (t: TooMany) => void,
): Map<string, PlannedTile> {
  const next = new Map<string, PlannedTile>();
  const planned: { w: Want; tiles: WantTile[]; shapes: readonly Shape[] }[] = [];
  const loads = new Map<string, LayerLoad & { keys: Set<string> }>();
  for (const w of wants) {
    const t = wantTiles(w, sources, settings);
    if (!t || t.count === 0) continue;
    if (t.count > settings.maxTilesPerWant) {
      tooMany({ want: w, bound: 'want', value: t.count, max: settings.maxTilesPerWant });
      continue;
    }
    const tiles = eachTile(w, t);
    planned.push({ w, tiles, shapes: t.shapes });
    if (w.layer === undefined) continue;
    let load = loads.get(w.layer);
    if (!load) loads.set(w.layer, (load = { tiles: 0, cost: 0, keys: new Set() }));
    const opened = sources.get(w.source)!.opened;
    for (const tile of tiles) {
      if (load.keys.has(tile.key)) continue;
      load.keys.add(tile.key);
      load.tiles++;
      if (opened.readCost) load.cost += costOf(opened, tile.region, t.shapes);
    }
  }
  for (const { w, tiles, shapes } of planned) {
    const load = w.layer === undefined ? undefined : loads.get(w.layer);
    if (load && load.tiles > settings.maxTilesPerLayer) {
      tooMany({ want: w, bound: 'layer', value: load.tiles, max: settings.maxTilesPerLayer });
      continue;
    }
    if (load && load.cost > settings.maxBytesPerLayer) {
      tooMany({ want: w, bound: 'cost', value: load.cost, max: settings.maxBytesPerLayer });
      continue;
    }
    const rank = Rank[w.priority ?? 'view'];
    for (const { key, region } of tiles) {
      const at = next.get(key);
      if (at) at.rank = Math.min(at.rank, rank);
      else next.set(key, { rank, region, shapes });
    }
  }
  return next;
}

/**
 * What a layer's wants cover together, and whether they are planned (each within `maxTilesPerWant`, together within
 * `maxTilesPerLayer` and `maxBytesPerLayer`): what a layer asks before it wants, to read coarser first.
 */
export function layerLoad(
  wants: readonly Want[],
  sources: StoreSources,
  settings: StoreSettings,
): LayerLoad & { fits: boolean } {
  const keys = new Set<string>();
  let cost = 0;
  let fits = true;
  for (const w of wants) {
    const t = wantTiles(w, sources, settings);
    if (!t || t.count === 0) continue;
    if (t.count > settings.maxTilesPerWant) {
      fits = false;
      continue;
    }
    const opened = sources.get(w.source)!.opened;
    for (const tile of eachTile(w, t)) {
      if (keys.has(tile.key)) continue;
      keys.add(tile.key);
      if (opened.readCost) cost += costOf(opened, tile.region, t.shapes);
    }
  }
  fits &&= keys.size <= settings.maxTilesPerLayer && cost <= settings.maxBytesPerLayer;
  return { tiles: keys.size, cost, fits };
}

/** What a want set asks of its store. */
export interface WantHost {
  _plan(wants: readonly Want[]): Map<string, PlannedTile>;
  _apply(): void;
  _removeSet(set: WantSet): void;
  _waitFor(check: () => boolean, timeoutMs?: number): Promise<boolean>;
  get(key: string): { readonly status: string } | undefined;
}

export class WantSet implements Wanter {
  tiles = new Map<string, PlannedTile>();
  /** As last asked: planned again when a source is read anew. */
  wants: readonly Want[] = [];
  private _disposed = false;

  constructor(private readonly _store: WantHost) {}

  get keys(): ReadonlySet<string> {
    return new Set(this.tiles.keys());
  }

  want(wants: readonly Want[]): void {
    if (this._disposed) return;
    this.wants = wants;
    this.tiles = this._store._plan(wants);
    this._store._apply();
  }

  whenSettled(timeoutMs?: number): Promise<boolean> {
    return this._store._waitFor(() => {
      for (const k of this.tiles.keys()) if (this._store.get(k)?.status === 'pending') return false;
      return true;
    }, timeoutMs);
  }

  dispose(): void {
    if (this._disposed) return;
    this._disposed = true;
    this.tiles = new Map();
    this.wants = [];
    this._store._removeSet(this);
  }
}
