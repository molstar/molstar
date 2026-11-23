/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PairsData, StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { LayoutContext, StoreReader } from './layer-type.js';
import { type PairsPrep, prepare } from './pairs-lay.js';
import { heldBins } from './pairs-held.js';
import { regionPairs, regionState } from './pairs-region.js';

// The pairs reader (_spec/README.md "Render", The pairs layer: the reader): one reader of a pairs node, held or read by
// region (one system or an assembly's blocks), giving the frame's pairs: bins at the pairs level (no wider than
// `binPx` device pixels), or the items themselves once the bin pairs with an end in view sum to at most its `gate`
// (kept to `hysteresis` × it). Held pairs past the gate are binned by the same level rule, so held and read pairs are
// one picture. The pairs layer reads through it with the pairs config (`pairs.binPx`, `pairs.maxArcs`,
// `pairs.itemsHysteresis`); an extension's layer type with its own numbers, its state kept by the reader.

/** What a read is asked: the bins' width in device pixels and the items gate; the pairs config's when absent. */
export interface PairsReadOptions {
  readonly binPx?: number;
  readonly gate?: number;
  readonly hysteresis?: number;
}

/**
 * A read's frame of pairs: its rows (ends on the layer's two systems, a pair found from its other end turned round;
 * over an assembly each end's unit, an index into `endUnits`), and for bins their level and counts.
 */
export interface PairsRead {
  /** The pairs node's address (a hit's element owner). */
  readonly owner: string;
  /** The layer's two systems: the ends' when no unit is given. */
  readonly aSystem: string;
  readonly bSystem: string;
  readonly data?: PairsData;
  /** 1 where the layer's data filter keeps a row. */
  readonly shown: Uint8Array;
  /** The bins' level in positions; absent: the rows are items. */
  readonly level?: number;
  /** A bin pair's count (absent: items, one each). */
  readonly weight?: Float32Array;
  /** Over an assembly: the sequences' systems the ends' units index. */
  readonly endUnits?: readonly string[];
  /** The rows the layer's data filter left out (held pairs). */
  readonly filtered: number;
}

/** What a reader reads: a pairs node's prep (opaque to the reader's owner), made once a resolved layer. */
export type PairsSource = PairsPrep & PairsRead;

/**
 * The frame's pairs of `prep` read for `owner` (whose state it keeps: a lay, a reader), or undefined while nothing is
 * read: a source by region (bins, then items under the gate), held pairs (binned past the gate, else themselves).
 */
export function readPairs(
  prep: PairsPrep,
  ctx: LayoutContext,
  data: StoreReader,
  owner: object,
  options?: PairsReadOptions,
): PairsPrep | undefined {
  if (prep.region && !prep.data) {
    const read = ctx.place
      ? regionPairs(prep, prep.region, ctx, data, owner, undefined, undefined, options)
      : undefined;
    return read?.data ? read : undefined;
  }
  if (!prep.data) return undefined;
  if (!prep.region && ctx.place) return heldBins(prep, ctx, data, options) ?? prep;
  return prep;
}

/**
 * A pairs reader for a layer type of its own (an extension's ribbons): `source` once a resolved layer (its `prepare`),
 * `read` once a layout with the type's numbers. One reader a lay: it keeps what it read between frames (a warm frame
 * builds nothing), the gate's state and the tiles pending.
 */
export class PairsReader {
  /** The pairs node a layer reads (its data, a source's region, its filter), made once a resolved layer. */
  static source(layer: ResolvedLayer, model: StaticModel): PairsSource {
    return prepare(layer, model) as PairsSource;
  }

  /** This frame's pairs, or undefined while nothing is read (a source's bins pending, a layout with no placement). */
  read(source: PairsSource, ctx: LayoutContext, data: StoreReader, options?: PairsReadOptions): PairsRead | undefined {
    this.last = readPairs(source, ctx, data, this, options) as PairsSource | undefined;
    return this.last;
  }

  /** The last read. */
  last: PairsRead | undefined = undefined;

  /** Whether the last read gave items (under the gate), not bins. */
  get items(): boolean {
    return this.last?.data !== undefined && this.last.level === undefined;
  }

  /** The bins' tiles not loaded as last read, `[start, end, …]` on the data's system (the pending hatch, EC-169). */
  get pending(): readonly number[] | undefined {
    return regionState(this).pending;
  }

  /** Whether the last read waits on a bins' tile it wanted (a region's or an assembly's): pending, not empty. */
  get waiting(): boolean {
    return regionState(this).waiting;
  }

  /** Whether a bins' tile of the last read failed (a source whose opening failed): failed, not pending. */
  get failed(): boolean {
    return regionState(this).failed;
  }
}
