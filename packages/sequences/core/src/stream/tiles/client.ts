/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `seqstar-tiles/1` (seqspec:sources#the-tiles-protocol): the client over a transport (the reference transport is HTTP
// GET under the source's uri). Binary payloads are read little-endian.

import { StreamNotServed, TILE, type Tile } from "@molstar/seq-core/format";

/**
 * The bytes of a resource under a source's uri (`tiles/0/0/0/0`). Rejects with `NotServed` when the source does not
 * serve it (HTTP 404 or 410), and with the failure otherwise: an outage, an abort.
 */
export type TileTransport = (path: string, signal?: AbortSignal) => Promise<Uint8Array>;

/** The source does not serve this resource: reported once, never retried. Only `stats` turns it into "none". */
export class NotServed extends StreamNotServed {
  constructor(message: string) {
    super(message);
    this.name = "NotServed";
  }
}

/** Binary payloads are little-endian: a view on a copy where the host is, a DataView read where it is not. */
const LITTLE_ENDIAN = new Uint8Array(Uint16Array.of(1).buffer)[0] === 1;
function floats(b: Uint8Array): Float64Array {
  if (LITTLE_ENDIAN) return new Float64Array(b.slice().buffer, 0, b.length >> 3);
  const v = new DataView(b.buffer, b.byteOffset, b.byteLength);
  return Float64Array.from({ length: b.length >> 3 }, (_, i) => v.getFloat64(8 * i, true));
}
function counts(b: Uint8Array): Uint32Array {
  if (LITTLE_ENDIAN) return new Uint32Array(b.slice().buffer);
  const v = new DataView(b.buffer, b.byteOffset, b.byteLength);
  return Uint32Array.from({ length: b.length >> 2 }, (_, i) => v.getUint32(4 * i, true));
}
const SAMPLE_BYTES = { tile: 2 * TILE * TILE, oneD: 16 * TILE, counts: 4 * TILE };
const decoder = new TextDecoder();

/** The client of one source: each call asks `transport` once. Rows 0-based, in leaf order. */
export interface TilesClient {
  readonly uri: string;
  tile(ly: number, lx: number, ty: number, tx: number, signal?: AbortSignal): Promise<Tile>;
  /** The labels of rows [a, b). */
  labels(a: number, b: number, signal?: AbortSignal): Promise<string[]>;
  /** Row `k` whole, the gapped string, and its label. */
  row(k: number, signal?: AbortSignal): Promise<{ label: string; row: string }>;
  /** A column statistic's 1-D tile, [min, max] pairs (NaN: none); null: not served. */
  stats(measure: string, lx: number, tx: number, signal?: AbortSignal): Promise<Float64Array | null>;
  /** An annotation's 1-D tiles: [min, max] pairs, or a count per sample. */
  values(lx: number, tx: number, signal?: AbortSignal): Promise<Float64Array>;
  counts(lx: number, tx: number, signal?: AbortSignal): Promise<Uint32Array>;
}

export function tileClient(uri: string, transport: TileTransport): TilesClient {
  const get = async (path: string, bytes: number, signal?: AbortSignal) => {
    const b = await transport(path, signal);
    if (bytes >= 0 && b.length !== bytes) throw new Error(`${uri}/${path}: ${b.length} bytes, not ${bytes}`);
    return b;
  };
  const json = async (path: string, signal?: AbortSignal) => JSON.parse(decoder.decode(await get(path, -1, signal)));
  /** A JSON answer checked as the binary ones are (GI-13): the protocol's shape, or an error naming the path. */
  const wrong = (path: string, what: string): never => {
    throw new Error(`${uri}/${path}: the answer is not ${what}`);
  };
  return {
    uri,
    async tile(ly, lx, ty, tx, signal): Promise<Tile> {
      const b = await get(`tiles/${ly}/${lx}/${ty}/${tx}`, SAMPLE_BYTES.tile, signal);
      return { cls: b.subarray(0, TILE * TILE), share: b.subarray(TILE * TILE) };
    },
    async labels(a, b, signal) {
      const path = `labels/${a}/${b}`;
      const v: unknown = await json(path, signal);
      const ok = Array.isArray(v) && v.length === b - a && v.every((x) => typeof x === "string");
      return ok ? (v as string[]) : wrong(path, `a list of ${b - a} labels`);
    },
    async row(k, signal) {
      const path = `rows/${k}`;
      const v = (await json(path, signal)) as { label?: unknown; row?: unknown } | null;
      const ok = v !== null && typeof v === "object" && typeof v.label === "string" && typeof v.row === "string";
      return ok ? (v as { label: string; row: string }) : wrong(path, "a row: {label, row}, both strings");
    },
    async stats(measure, lx, tx, signal) {
      return get(`stats/${measure}/${lx}/${tx}`, SAMPLE_BYTES.oneD, signal).then(floats, (e) => {
        if (e instanceof NotServed) return null;
        throw e;
      });
    },
    values: async (lx, tx, signal) => floats(await get(`values/${lx}/${tx}`, SAMPLE_BYTES.oneD, signal)),
    async counts(lx, tx, signal) {
      return counts(await get(`counts/${lx}/${tx}`, SAMPLE_BYTES.counts, signal));
    },
  };
}

/** The reference transport: HTTP GET under `uri` (global `fetch`). */
export function fetchTransport(uri: string): TileTransport {
  return async (path, signal) => {
    const r = await fetch(`${uri.replace(/\/$/, "")}/${path}`, { signal });
    if (r.status === 404 || r.status === 410) throw new NotServed(`${r.status} ${r.statusText}`);
    if (!r.ok) throw new Error(`${r.status} ${r.statusText}`);
    return new Uint8Array(await r.arrayBuffer());
  };
}
