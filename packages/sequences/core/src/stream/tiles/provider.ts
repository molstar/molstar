/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `seqstar-tiles/1` on the stream port (seqspec:sources#the-tiles-protocol, "The provider"): one provider per source,
// built from the shape its readers declare (`StreamSource.declared`): a manifest, an alignment's tiles and rows, an
// annotation's bins.

import {
  type DataPatch,
  type Manifest,
  StreamNotServed,
  type StreamProvider,
  type StreamRequest,
  type StreamSource,
  TILE,
  Tiles,
} from "@molstar/seq-core/format";
import { fetchTransport, type TilesClient, type TileTransport, tileClient } from "./client.js";

export const PROTOCOL = "seqstar-tiles/1";
/** An alignment's levels: every [2^ly, 2^lx], level `ly·16 + lx`. */
const GRID: [number, number][] = Array.from({ length: 256 }, (_, k) => [2 ** (k >> 4), 2 ** (k & 15)]);

/**
 * A host's `loaders.stream` for `seqstar-tiles/1` sources: `open` gives the transport for a source's uri (null: not
 * served here), by default HTTP GET under it. Any other source gets null.
 */
export function tilesStream(open: (uri: string) => TileTransport | null = fetchTransport) {
  return (source: StreamSource): StreamProvider | null => {
    const transport = source.protocol === PROTOCOL && source.uri != null ? open(source.uri) : null;
    return transport && tilesProvider(source, tileClient(source.uri as string, transport));
  };
}

/**
 * A source's provider over its client. Its readers' params are in `source.declared`: an alignment's `row_count` and
 * `column_count`, an annotation's `length` and `target` (its one coordinate system's key), and `measure` for a
 * column statistic (a `values` source the alignment's `column_stats` implies). A manifest request for a statistic the
 * server does not serve rejects with `StreamNotServed`.
 */
export function tilesProvider(source: StreamSource, client: TilesClient): StreamProvider {
  const kind = source.data_kind;
  const d = source.declared;
  const measure = typeof d.measure === "string" ? d.measure : null;
  const declared = (k: string) => {
    if (typeof d[k] !== "number") throw new Error(`${client.uri}: a ${kind} source declares ${k}`);
    return d[k] as number;
  };
  const read = async (lx: number, tx: number, signal: AbortSignal) => {
    if (kind === "items") return client.counts(lx, tx, signal);
    const t = await (measure ? client.stats(measure, lx, tx, signal) : client.values(lx, tx, signal));
    if (!t) throw new StreamNotServed(`${client.uri} serves no stats/${measure}`);
    return t;
  };
  async function manifest(signal: AbortSignal): Promise<Manifest> {
    const system = (extent: number) => [
      { key: String(d.target ?? source.ref), extent: [0, extent] as [number, number], numbering: null, label: null },
    ];
    if (kind === "alignment") {
      const [rows, columns] = [declared("row_count"), declared("column_count")];
      return { systems: system(columns), levels: GRID, rows, columns, labels: true };
    }
    if (kind !== "values" && kind !== "items") throw new StreamNotServed(`${PROTOCOL} serves no ${kind} sources`);
    const length = declared("length");
    const top = Tiles.levelFor(length);
    const levels = Array.from({ length: top + 1 }, (_, k) => 2 ** k);
    if (kind === "items") return { systems: system(length), levels, vocabulary: {} };
    // The domain: the range of the coarsest tile (one sample a 2^top positions).
    const t = await read(top, 0, signal);
    let [lo, hi] = [Infinity, -Infinity];
    for (let i = 0; i < t.length; i += 2) {
      if (t[i] < lo) lo = t[i];
      if (t[i + 1] > hi) hi = t[i + 1];
    }
    if (!(lo <= hi)) throw new StreamNotServed(`${client.uri}: ${measure ?? "values"} has no value`);
    return { systems: system(length), levels, domain: [lo, hi] };
  }
  /** A region's bins: of `request.level`'s width, else (the data) of one position, its tiles asked at once, stitched. */
  async function bins(r: StreamRequest, signal: AbortSignal): Promise<Record<string, ArrayLike<number>>> {
    const lx = Math.max(0, r.level);
    const w = 2 ** lx;
    const n = Math.ceil((r.stop - r.start) / w);
    const first = Math.floor(r.start / w);
    const t0 = Math.floor(first / TILE);
    const asked: Promise<Float64Array | Uint32Array>[] = [];
    for (let tx = t0; tx <= Math.floor((first + n - 1) / TILE); tx++) asked.push(read(lx, tx, signal));
    const tiles = await Promise.all(asked);
    if (kind === "items") {
      const count = new Uint32Array(n);
      for (let j = 0; j < n; j++) count[j] = tiles[Math.floor((first + j) / TILE) - t0][(first + j) % TILE];
      return { count };
    }
    const min = new Float64Array(n);
    const max = new Float64Array(n);
    const mean = new Float64Array(n);
    for (let j = 0; j < n; j++) {
      const t = tiles[Math.floor((first + j) / TILE) - t0];
      const s = (first + j) % TILE;
      min[j] = t[2 * s];
      max[j] = t[2 * s + 1];
      mean[j] = (min[j] + max[j]) / 2;
    }
    return { min, max, mean };
  }
  return {
    async fetch(r: StreamRequest, signal: AbortSignal): Promise<DataPatch> {
      const at = { source: r.source, system: r.system, start: r.start, stop: r.stop, level: r.level };
      const rows = { row_start: r.row_start, row_stop: r.row_stop };
      if (r.system === null) return { ...at, ...rows, kind: "manifest", manifest: await manifest(signal) };
      if (kind !== "alignment") return { ...at, ...rows, kind: "bins", bins: await bins(r, signal) };
      if (r.labels) {
        const labels = await client.labels(r.row_start, r.row_stop, signal);
        return { ...at, ...rows, kind: "rows", labels, rows: null };
      }
      if (r.level < 0) {
        const got = await Promise.all(
          Array.from({ length: r.row_stop - r.row_start }, (_, j) => client.row(r.row_start + j, signal)),
        );
        return { ...at, ...rows, kind: "rows", labels: got.map((x) => x.label), rows: got.map((x) => x.row) };
      }
      const [ly, lx] = [r.level >> 4, r.level & 15];
      const t = await client.tile(ly, lx, Tiles.tileOf(r.row_start, ly), Tiles.tileOf(r.start, lx), signal);
      return { ...at, ...rows, kind: "bins", bins: { class: t.cls, share: t.share } };
    },
  };
}
