/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What BigBed and tabix VCF share (README, "Items, levels and the zoom gate"): an `items` provider over a reader of
// the items that meet a region. Levels are counts per bin, computed from the items or served by the source without
// reading them (ST4b: a BigBed's zoom records, a tabix index's estimate); a data request holding more than `maxItems`
// (estimated first where the source can, so it is not read) is answered with the counts of the finest level that fits
// (the zoom gate, seqspec:sources#streamed-sources).

import {
  type DataPatch,
  type Manifest,
  StreamNotServed,
  type StreamProvider,
  type StreamRequest,
  TILE,
} from "@molstar/seq-core/format";

/** Positions per bin of the levels an items source lists (a level's region is `TILE` of them: 64 kb, 1 Mb, 16 Mb). */
export const ITEM_LEVELS = [256, 4096, 65536];
export const MAX_ITEMS = 2000;

export interface System {
  readonly key: string;
  readonly extent: number;
}
/** The items meeting a region, whole: `table` has a row per drawn piece, `start` and `stop` one entry per item. */
export interface ItemsRead {
  readonly table: Record<string, ArrayLike<unknown>>;
  readonly start: readonly number[];
  readonly stop: readonly number[];
}
export interface ItemsSource {
  systems(signal: AbortSignal): Promise<System[]>;
  read(system: string, start: number, stop: number, signal: AbortSignal): Promise<ItemsRead>;
  vocabulary: Record<string, string[]>;
  /** Positions per bin of its levels, finest first, read with its manifest (default `ITEM_LEVELS`; ST4b). */
  levels?(signal: AbortSignal): Promise<readonly number[]>;
  /**
   * Level `level`'s bins over [start, stop) without reading the items (a BigBed's zoom records, a tabix index's
   * estimate; ST4b); null: they are counted from the items (a fine level, a region cheap to read).
   */
  bins?(system: string, start: number, stop: number, level: number, signal: AbortSignal): Promise<Float64Array | null>;
  /** Its `bins` are estimates (an index's), not counts: the patch says so (`estimate`). */
  readonly estimated?: boolean;
  /** What reading [start, stop) would cost, from an index: its items and bytes (ST4b); null: unknown, it is read. */
  estimate?(system: string, start: number, stop: number, signal: AbortSignal): Promise<Estimate | null>;
}
export interface Estimate {
  readonly items: number;
  readonly bytes: number;
}

/**
 * An `items` provider over `source`. A level's bins come from the source where it serves them without reading the
 * items, else from the items. A data request is read unless the source's estimate says it holds more than `maxItems`
 * items or more bytes than the request's `max_bytes`, or once read it holds more than `maxItems`: then it is answered
 * with the bins of the finest level that fits a tile (the zoom gate); a `force` request is read whole.
 */
export function itemsProvider(source: ItemsSource, maxItems = MAX_ITEMS): StreamProvider {
  let levels: readonly number[] | null = null;
  const levelsOf = async (signal: AbortSignal) => {
    levels ??= source.levels ? await source.levels(signal) : ITEM_LEVELS;
    return levels;
  };
  return {
    async fetch(r: StreamRequest, signal: AbortSignal): Promise<DataPatch> {
      const at = {
        source: r.source,
        system: r.system,
        start: r.start,
        stop: r.stop,
        row_start: 0,
        row_stop: 0,
        level: r.level,
      };
      const L = await levelsOf(signal);
      if (r.system === null) {
        const systems = (await source.systems(signal)).map((s) => ({
          key: s.key,
          extent: [0, s.extent] as [number, number],
          numbering: null,
          label: s.key,
        }));
        const manifest: Manifest = { systems, levels: [...L], vocabulary: source.vocabulary };
        return { ...at, kind: "manifest", manifest };
      }
      if (r.level >= L.length) throw new StreamNotServed(`no level ${r.level}: this source lists ${L.length}`);
      if (r.start < 0 || r.stop < r.start) throw new Error(`bad region [${r.start}, ${r.stop})`);
      const binned = (level: number, count: ArrayLike<number>, estimate = false): DataPatch => ({
        ...at,
        level,
        kind: "bins",
        bins: { count },
        ...(estimate ? { estimate } : {}),
      });
      // The finest level whose bins for this region fit in a tile: what a gated data request is answered with.
      const fit = L.findIndex((w) => (r.stop - r.start) / w <= TILE);
      const gate = r.level >= 0 ? r.level : fit < 0 ? L.length - 1 : fit;
      const served = (level: number) =>
        source.bins ? source.bins(r.system as string, r.start, r.stop, level, signal) : Promise.resolve(null);
      if (r.level >= 0) {
        const pre = await served(r.level);
        if (pre) return binned(r.level, pre, source.estimated === true);
      } else if (!r.force && source.estimate) {
        const e = await source.estimate(r.system, r.start, r.stop, signal);
        const over = e !== null && (e.items > maxItems || (r.max_bytes != null && e.bytes > r.max_bytes));
        const pre = over ? await served(gate) : null;
        if (pre) return binned(gate, pre, source.estimated === true);
      }
      const got = await source.read(r.system, r.start, r.stop, signal);
      if (r.level < 0 && (r.force || got.start.length <= maxItems)) return { ...at, kind: "items", items: got.table };
      const w = L[gate];
      const count = new Uint32Array(Math.ceil((r.stop - r.start) / w));
      for (let i = 0; i < got.start.length; i++) {
        const last = Math.min(count.length - 1, Math.floor((got.stop[i] - 1 - r.start) / w));
        for (let j = Math.max(0, Math.floor((got.start[i] - r.start) / w)); j <= last; j++) count[j]++;
      }
      return binned(gate, count);
    },
  };
}
