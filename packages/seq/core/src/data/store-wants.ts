/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { dataKey, fitLevel, tileRegion, tileSpan } from './keys.js';
import type { Shape } from './protocol.js';
import type { StoreSources } from './store-sources.js';
import { Rank, type StoreSettings, type TileRegion, type Want, type Wanter } from './store-types.js';
import { sourceNameOf } from './system-names.js';

// The data store's want sets and their planning (_spec/README.md "Data", Store): a consumer's wants made into the
// tiles they cover, each at its best priority; a want past `maxTilesPerWant` is not planned (said by the caller).

/** A tile a want set needs, at its best rank. */
export interface PlannedTile {
  rank: number;
  readonly region: TileRegion;
  readonly shapes: readonly Shape[];
}

/** The tiles of a want: its system's coverage name, the level read, the tiles' span; undefined for an unknown one. */
export function wantTiles(w: Want, sources: StoreSources, settings: StoreSettings) {
  const s = sources.get(w.source);
  const system = s && sourceNameOf(s.names, w.system);
  if (!s || system === undefined) return undefined;
  const cov = s.opened.coverage;
  const level = cov.levels.includes(w.level) ? w.level : fitLevel(cov.levels, w.level);
  const length = cov.systems.find((x) => x.name === system)?.length;
  const end = length === undefined ? w.end : Math.min(w.end, length);
  const [first, last] = tileSpan(w.start, end, level, settings.tileBins);
  const shapes = cov.shapes.filter((x) => x.system === system).map((x) => x.shape);
  return { system, level, length, first, last, shapes };
}

/** A consumer's wants as tiles by key; `tooMany` is told of each want past `maxTilesPerWant` (not planned). */
export function planWants(
  wants: readonly Want[],
  sources: StoreSources,
  settings: StoreSettings,
  tooMany: (w: Want, tiles: number) => void,
): Map<string, PlannedTile> {
  const next = new Map<string, PlannedTile>();
  for (const w of wants) {
    const t = wantTiles(w, sources, settings);
    if (!t || t.last < t.first) continue;
    if (t.last - t.first + 1 > settings.maxTilesPerWant) {
      tooMany(w, t.last - t.first + 1);
      continue;
    }
    const rank = Rank[w.priority ?? 'view'];
    for (let tile = t.first; tile <= t.last; tile++) {
      const key = dataKey(w.source, t.system, t.level, tile);
      const at = next.get(key);
      if (at) at.rank = Math.min(at.rank, rank);
      else {
        const { start, end } = tileRegion(tile, t.level, settings.tileBins, t.length);
        const region = { source: w.source, system: t.system, level: t.level, tile, start, end };
        next.set(key, { rank, region, shapes: t.shapes });
      }
    }
  }
  return next;
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
