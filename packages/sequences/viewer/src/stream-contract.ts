/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the stream manager promises its readers (seqspec:sources#the-viewers-side; seqspec:viewer#extension-points): a
// region's state, a layer's description, what drawing reads (`StreamView`) and the host's manager (`StreamManager`).
// Types and constants only; `streams.ts` re-exports them.

import type { StreamedData, Tile } from '@molstar/seq-core/format';
import type { ComposedPath } from '@molstar/seq-core/model';
import type { Features } from './prepare.js';
import type { Live, SeqState } from './state.js';

/** A region's state; one not held is pending (never drawn or hit as empty). */
export const RegionState = { pending: 0, loaded: 1, empty: 2, failed: 3 } as const;

/** An alignment in memory: its held rows' codes, and per lane its slot in them (-1: not held). */
export interface LocalCells {
  readonly codes: Uint8Array;
  readonly slot: Int32Array;
}
/** A layer's data held by region (`renderer.stream(s)`, seqspec:viewer#extension-points): a pure description. */
export interface StreamLayer {
  /** The data's identity, its content included: a layer's regions survive a re-prepare with the same key. */
  readonly key: string;
  /** An alignment's lanes and columns (0 for a 1-D layer). */
  readonly rows: number;
  readonly columns: number;
  /** A streamed source: its provider, manifest and labels (`svs:stream`). */
  readonly remote: StreamedData | null;
  readonly local: LocalCells | null;
  /** 1-D: per block of its section, the manifest system it shows (-1: none); null: an alignment's tiles. */
  readonly systems: Int32Array | null;
  /** Items: a field whose equal values fold the items into one lane (a `ranges` layer's `collapse`); null or absent: none. */
  readonly collapse?: string | null;
  /**
   * 1-D, per block: the path from its system to the block's (ST7: a genome's data on a protein's block); null or
   * absent: every block shows its system directly.
   */
  readonly paths?: readonly (ComposedPath | null)[] | null;
}
export interface StreamStats {
  requests: number;
  hits: number;
  inFlight: number;
  cancelled: number;
  /** Arrivals of coarse (and overview) regions, and of the view's level. */
  coarse: number;
  fine: number;
  /** Held by the cache. */
  bytes: number;
  errors: number;
  /** Local tiles finished. */
  computed: number;
}
/** `Stream.LettersPerPx`'s default: a letters layer asks its data only where a device px holds at most this many. */
export const LETTERS_PER_PX = 16;
/**
 * One region's items (`StreamView.items`): a feature table built at arrival (over the section's blocks that show its
 * system, lanes the layer's by item, fragments of one `id` joined), or a level's `count` bins.
 */
export interface StreamItems {
  /** Its items; null: bins, or nothing held. */
  readonly features: Features | null;
  /** Per feature its item's start: an item is drawn by the first region in view that holds it. */
  readonly first: Int32Array;
  /**
   * Per feature the start of what it was drawn from, for that rule, when it differs from `first` (a fold's group
   * begins where its first member does, but a member is held by every region it overlaps): absent: `first`.
   */
  readonly own?: Int32Array;
  /** A level's `count` per bin from the region's start (density); null: items. */
  readonly bins: Float64Array | null;
  /** Positions per bin. */
  readonly width: number;
  /** The patch's columns (`row` indexes them): fields for labels and field themes; null: bins. */
  readonly table: Readonly<Record<string, ArrayLike<unknown>>> | null;
  /** Per feature its item's number in the layer (the hover's `part`); null: bins. */
  readonly item: Int32Array | null;
  /** Per feature its item's stop, as far as the region's patch holds it (`first` its start): the item's range. */
  readonly ends: Int32Array | null;
  /** Bins a reader should be told about (ST4b): 1, an index's estimate (zoom in); 2, a gated data region (zoom in, or
   * load it anyway); 0: none. */
  readonly hint: number;
}
/** What drawing reads of a layer (`DrawView.stream`): its regions by level and tile. Allocation-free. */
export interface StreamView {
  readonly source: StreamLayer;
  /** Bumped when a region of the layer arrives, fails or is evicted. */
  readonly version: number;
  /**
   * `Stream.SampleBudget` as the manager read it, so draws pick the manager's levels; 0 for null: a sample per device
   * pixel at most, which the levels' rounding gives alone (a small integer: an infinity passed to a call is boxed).
   */
  readonly budget: number;
  /** `Stream.LettersPerPx` as the manager read it (the letters gate; none: `LETTERS_PER_PX`). */
  readonly gate?: number;
  /**
   * 2-D: per level `ly * 16 + lx` the level its source serves in its place (`listedLevels`; -1: it serves none); null or
   * absent: every level (a local stream).
   */
  readonly served?: Int16Array | null;
  tile(ly: number, lx: number, ty: number, tx: number): Tile | null;
  state(ly: number, lx: number, ty: number, tx: number): number;
  /** A row's label once its source served it. */
  label(row: number): string | null;
  /**
   * 1-D `values`: block b's indices [i0, i1) at `out[0]` positions a device px (read first; a double crosses no call)
   * into `out` the minimum, maximum, sum of means and count of the bins whose start lies in it (else the bin holding
   * i0), from the view's level, else a coarser one loaded (never blank). False: none loaded.
   */
  values?(b: number, i0: number, i1: number, out: Float64Array): boolean;
  /**
   * 1-D `letters`: block b's indices [i0, i1) into `out[0, i1 − i0)` (at most `out.length`), the data's upper-case codes,
   * 0 where no region holding them is loaded; how many are held (0: draw nothing).
   */
  letters?(b: number, i0: number, i1: number, out: Uint8Array): number;
  /**
   * 1-D `items`: the k-th region (from 0, left to right) over block b's indices [i0, i1) at `out[0]` positions a device
   * px (the view's level; -1: the data); null past the last. Into `out[1..5]`: its [start, stop) in block indices, the position below
   * which its items are an earlier region's (an item crossing regions is held in each and drawn once, by the first),
   * its `RegionState` and the layer's lanes; where `out` has room, `out[6]` its start in its system's positions (where
   * its bins and its items' `first` and `ends` are). A region not loaded gives no features and no bins.
   */
  items?(b: number, i0: number, i1: number, k: number, out: Int32Array): StreamItems | null;
  /**
   * 1-D: positions [s0, s1) of block b's system as the hull of the block's indices they reach into `out[0, 2)` (through
   * its path; else themselves); false: none.
   */
  toBlock?(b: number, s0: number, s1: number, out: Int32Array): boolean;
  /**
   * 1-D `items`, a draw's, once a frame before its `items` reads: the lanes that fit the band. A layer with a `collapse`
   * field gives its items folded (one lane a group, the parts of a group's members) from then on while its lanes would
   * not fit and its track is not expanded; else its own.
   */
  fold?(fit: number): void;
  /** 1-D: the `RegionState` of block b's index i at the view's level (past the letters gate, out of the source: empty). */
  stateAt?(b: number, i: number): number;
  /**
   * 1-D, the overlay's: the k-th region at the view's level over block b's indices [i0, i1), `out[0]` the frame's time
   * (ms): into `out[1..4]` its [start, stop) in block indices, its `RegionState` (values: loaded where a coarser level is
   * drawn in its place) and its age, the ms since the frame that first drew it (its coarser level's, values); false past
   * the last, or none needed there (out of view, past the letters gate).
   */
  region?(b: number, i0: number, i1: number, k: number, out: Float64Array): boolean;
  /** 2-D, the overlay's: a tile's `RegionState`, its age into `out[1]` as `region`'s (`out[0]` the frame's time). */
  tileAt?(ly: number, lx: number, ty: number, tx: number, out: Float64Array): number;
  /**
   * The source columns [c0, c1) of the regions that changed after version `since`, pairs into `out`: their count; -1
   * when not known (too many, a labels chunk, or `out` too small): everything.
   */
  changed?(since: number, out: Int32Array): number;
}
export interface StreamManager {
  get(layer: string): StreamView | null;
  /**
   * RegionState of the layer's region holding index `index` of block `block`, at the level the view needs; a rows
   * layer's at `lane` (the layer's own lane, from 0; default the first in view).
   */
  state(layer: string, block: number, index: number, lane?: number): number;
  /** Every region in view loaded or failed at the level it needs (an inactive or disposed manager: now). */
  settled(): Promise<void>;
  /** Ends the cooldown of the layer's failed regions (not `NotServed` ones) and requests them; how many. */
  retry(layer: string): number;
  /**
   * The layer's gated data region at block `block`'s index `index` asked again whole (a "Load anyway", the zoom gate
   * aside; ST4b); false: none there.
   */
  force(layer: string, block: number, index: number): boolean;
  /** True while a region of the layer failed and is not `NotServed`: what `retry` can end (the track menu's "Retry"). */
  hasFailed(layer: string): boolean;
  /**
   * A lazy molecule's letters (its stream `data`, manifest system `s`) at `n` positions `at` into `out` (0: not held,
   * or `at` -1); how many are held; -1: a region holding one failed (asked again after its backoff, or by `retry`),
   * -2: one is not served. Those not held are asked at the next pass after the view's regions, without blocking, each
   * read replacing the last one's (ST7: a hovered residue's codon). A view's letters layer of the stream is read first.
   * Allocation-free once its regions are held.
   */
  lettersAt?(data: StreamedData, s: number, at: Int32Array, n: number, out: Uint8Array): number;
  readonly version: SeqState<number>;
  readonly stats: Live<StreamStats>;
}
