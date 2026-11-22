/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue, Node } from '@molstar/core/repr/tree';
import { DefaultGridTileBins, DefaultTileBins } from './keys.js';
import type { FailureKind, Tables } from './protocol.js';
import { DefaultSummarySettings, type SummarySettings } from './summaries.js';
import { DefaultAnswerBounds, DefaultMaxBlockValues, DefaultMaxColumnLength, DefaultMaxRows } from './tables.js';

// The data store's types and settings (_spec/README.md "Data", Store): what a consumer wants, what an entry holds, the
// store's clock and its settings.

/** Where a wanted region stands: in the view, its margin, or prefetched (EC-172's order). */
export type Priority = 'view' | 'margin' | 'prefetch';

export const Rank: Record<Priority, number> = { view: 0, margin: 1, prefetch: 2 };
export const Unwanted = 3;

/** A region a visible layer needs: positions `[start, end)` of a system of an attached source, at a level. */
export interface Want {
  readonly source: string;
  readonly system: string;
  readonly start: number;
  readonly end: number;
  /** A bin width; an unlisted one falls back to the coarsest listed no wider (EC-170). */
  readonly level: number;
  /** `view` by default. */
  readonly priority?: Priority;
  /** A grid's second axis: the want's tiles are its x tiles × its y tiles (a symmetric grid's folded). */
  readonly y?: { readonly system: string; readonly start: number; readonly end: number };
  /** A sided pairs source's side read. */
  readonly side?: 'a' | 'b';
  /** Items tiles (`RegionRequest.items`): the level is only the tile's width. */
  readonly items?: boolean;
  /**
   * The layer the want is for: a layer's wants together are bounded (`maxTilesPerLayer`, `maxBytesPerLayer`); none, the
   * want alone.
   */
  readonly layer?: string;
}

/** The tile an entry holds: its parts and its positions. */
export interface TileRegion {
  readonly source: string;
  /** The coverage's name of the system (an alias resolved). */
  readonly system: string;
  readonly level: number;
  readonly tile: number;
  readonly start: number;
  readonly end: number;
  /** A 2D tile's second axis: its system, tile and positions (a symmetric grid's served tile, never its mirror). */
  readonly y?: { readonly system: string; readonly tile: number; readonly start: number; readonly end: number };
  readonly side?: 'a' | 'b';
  readonly items?: boolean;
}

/**
 * A data entry: pending is never drawn or hit as empty (EC-169); failed carries its issue; evicted is wanted but
 * dropped for the byte bound, not asked again while it stays wanted.
 */
export type DataEntry =
  | { readonly status: 'pending'; readonly version: number; readonly region: TileRegion }
  | {
      readonly status: 'loaded';
      readonly version: number;
      readonly region: TileRegion;
      readonly tables: Tables;
      /** The node keys (`<source>#<key>`) this region carried. */
      readonly nodes: readonly string[];
      readonly bytes: number;
    }
  | { readonly status: 'empty'; readonly version: number; readonly region: TileRegion }
  | {
      readonly status: 'failed';
      readonly version: number;
      readonly region: TileRegion;
      readonly kind: FailureKind;
      readonly issue: Issue;
    }
  | { readonly status: 'evicted'; readonly version: number; readonly region: TileRegion };

export type EntryStatus = DataEntry['status'];

/** A node entry: a dynamic node of a source, held while something holds it. */
export interface NodeEntry {
  readonly source: string;
  readonly key: string;
  readonly node: Node;
  readonly version: number;
}

/** The store's time: real timers by default, a fake one in tests. */
export interface StoreClock {
  now(): number;
  setTimeout(fn: () => void, ms: number): unknown;
  clearTimeout(handle: unknown): void;
}

export const RealClock: StoreClock = {
  now: () => performance.now(),
  setTimeout: (fn, ms) => setTimeout(fn, ms),
  clearTimeout: (h) => clearTimeout(h as ReturnType<typeof setTimeout>),
};

export interface StoreSettings {
  /** Bins per tile. */
  readonly tileBins: number;
  /** The data entries' byte bound. */
  readonly maxBytes: number;
  /** Requests in flight at most (EC-172). */
  readonly maxInFlight: number;
  /** Bins a side of a 2D (grid) tile: 256. */
  readonly gridTileBins: number;
  /** Tiles one want may cover (a 2D want's x tiles × y tiles); more is an issue and the want is not requested. */
  readonly maxTilesPerWant: number;
  /** Tiles a layer's wants may cover together: 1,024; past it the layer reads coarser. */
  readonly maxTilesPerLayer: number;
  /** A layer's tiles' read cost together (`OpenedSource.readCost`), bytes: 16 MiB; past it the layer reads coarser. */
  readonly maxBytesPerLayer: number;
  /** The bytes one answer may hold (`answerBytes`): 16 MB, a quarter of the cache; more is `too-large`. */
  readonly maxTileBytes: number;
  /** Rows a list or density table may hold in one answer. */
  readonly maxRows: number;
  /** The longest column an answer may hold: 2^31 − 1. */
  readonly maxColumnLength: number;
  /** The most values (rows times columns) one table of an answer may hold: 2^24. */
  readonly maxBlockValues: number;
  /** Dynamic nodes in one answer: 10,000. */
  readonly maxNodes: number;
  /** What the nodes of one answer weigh together (the spec's size estimate): 16 MB. */
  readonly maxNodeBytes: number;
  /** Issues kept; then one `seq.store.limit_exceeded`. */
  readonly maxIssues: number;
  /** How long a request that left the wanted set stays in flight before it is cancelled (EC-173). */
  readonly graceMs: number;
  /**
   * A `failed` answer: tried `attempts` times with backoff `baseMs · factor^n`, then a cooldown (a `malformed` one goes
   * to the cooldown at once).
   */
  readonly retry: {
    readonly attempts: number;
    readonly baseMs: number;
    readonly factor: number;
    readonly cooldownMs: number;
  };
  /** The value pyramid given a source that serves level 1 only (`withSummaries`); `false`: none. */
  readonly summaries: SummarySettings | false;
}

export const DefaultStoreSettings: StoreSettings = {
  tileBins: DefaultTileBins,
  maxBytes: 64 * 1024 * 1024,
  maxInFlight: 6,
  gridTileBins: DefaultGridTileBins,
  maxTilesPerWant: 256,
  maxTilesPerLayer: 1024,
  maxBytesPerLayer: 16 * 1024 * 1024,
  maxTileBytes: 16 * 1024 * 1024,
  maxRows: DefaultMaxRows,
  maxColumnLength: DefaultMaxColumnLength,
  maxBlockValues: DefaultMaxBlockValues,
  maxNodes: DefaultAnswerBounds.maxNodes,
  maxNodeBytes: DefaultAnswerBounds.maxNodeBytes,
  maxIssues: 100,
  graceMs: 300,
  retry: { attempts: 3, baseMs: 500, factor: 2, cooldownMs: 30_000 },
  summaries: DefaultSummarySettings,
};

/** A consumer's want set (the live view, an export): merged with the others by best priority. */
export interface Wanter {
  /** The regions this consumer needs now, all of them. */
  want(wants: readonly Want[]): void;
  /** `true` when none of this consumer's tiles is pending; `false` after `timeoutMs` (if given) without (EC-125). */
  whenSettled(timeoutMs?: number): Promise<boolean>;
  /** The keys of the entries this consumer wants. */
  readonly keys: ReadonlySet<string>;
  /** Withdraws its wants. */
  dispose(): void;
}
