/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The reference `seqstar-tiles/1` server (seqspec:sources#the-tiles-protocol) as an in-process transport: rows
// generated from their index unless given, tiles aggregated as the spec says (the most frequent residue class, ties to
// the lower code, `-` only where a sample has no residue, `.` read as `-`; its share of all the sample's cells × 255),
// 1-D tiles of [min, max] float64 pairs or uint32 counts. Every request is recorded in `calls`. It holds nothing: a
// procedural cell is a hash of its place, a tile is made when asked and a sample past `EXACT` cells reads a strided
// subset of them, so an alignment of a million rows costs what a tile does (ST2; the first version cached every row
// string, about 4 GB at a million rows).

import { TILE } from '@molstar/seq-core/format';
import { NotServed, type TileTransport } from './client.js';

export interface TileServerOptions {
  rowCount: number;
  columnCount: number;
  /** Row k's gapped string (upper case, `-` a gap); default: procedural. */
  row?: (k: number) => string;
  label?: (k: number) => string;
  /** Column statistics served (`stats/<measure>/`): a value per column. */
  stats?: Record<string, (column: number) => number>;
  /** A progressive annotation of `length` positions: `values/` and `counts/`. */
  length?: number;
  values?: (i: number) => number;
  counts?: (i: number) => number;
}

const LETTERS = 'ACDEFGHIKLMNPQRSTVWY';
/** The code of a procedural cell: a letter from a hash of (row, column), a gap (45) at about one cell in six. */
export function proceduralCode(k: number, c: number): number {
  const h = Math.imul((k + 1) * 73856093 + c * 19349663, 2654435761) >>> 0;
  return h % 6 === 0 ? 45 : LETTERS.charCodeAt((h >>> 8) % LETTERS.length);
}
/** A deterministic gapped row (`proceduralCode`'s cells). */
export function proceduralRow(k: number, columns: number): string {
  let out = '';
  for (let c = 0; c < columns; c++) out += String.fromCharCode(proceduralCode(k, c));
  return out;
}
/** The residues of a procedural row, counted from its cells (no string of a million letters to make). */
function residues(k: number, columns: number): number {
  let n = 0;
  for (let c = 0; c < columns; c++) if (proceduralCode(k, c) !== 45) n++;
  return n;
}

/** Cells a sample reads exactly; a sample of more reads about this many, strided over its rows and columns. */
const EXACT = 1024;

export function tileServer(o: TileServerOptions): { transport: TileTransport; calls: string[] } {
  const calls: string[] = [];
  const text = (k: number) => o.row?.(k) ?? proceduralRow(k, o.columnCount);
  const label =
    o.label ?? ((k: number) => `SEQ${k}/1-${o.row ? o.row(k).replace(/-/g, '').length : residues(k, o.columnCount)}`);
  const json = (v: unknown) => new TextEncoder().encode(JSON.stringify(v));
  const tile = (ly: number, lx: number, ty: number, tx: number) => {
    const out = new Uint8Array(2 * TILE * TILE);
    const count = new Uint32Array(128);
    // Given rows are asked once a tile call (a sample's rows are shared by the samples across), and let go with it;
    // read once a sampled row, not once a cell.
    const given = new Map<number, string>();
    const rowOf = (r: number) => {
      let s = given.get(r);
      if (s === undefined) {
        s = (o.row as (k: number) => string)(r);
        given.set(r, s);
      }
      return s;
    };
    for (let sy = 0; sy < TILE; sy++) {
      for (let sx = 0; sx < TILE; sx++) {
        const [r0, c0] = [(ty * TILE + sy) * 2 ** ly, (tx * TILE + sx) * 2 ** lx];
        const [r1, c1] = [Math.min(r0 + 2 ** ly, o.rowCount), Math.min(c0 + 2 ** lx, o.columnCount)];
        if (r0 >= r1 || c0 >= c1) continue;
        // About EXACT cells at most, by a stride over rows and columns: a million rows are sampled, not read.
        const extra = ((r1 - r0) * (c1 - c0)) / EXACT;
        const dr = extra > 1 ? Math.min(r1 - r0, Math.ceil(Math.sqrt(extra))) : 1;
        const dc = extra > 1 ? Math.ceil(extra / dr) : 1;
        count.fill(0);
        let n = 0;
        // `.` is a gap as `-` is (the parser normalizes it, so the in-memory reduction never sees it).
        for (let r = r0; r < r1; r += dr) {
          const s = o.row ? rowOf(r) : null;
          for (let c = c0; c < c1; c += dc, n++) {
            const k = s !== null ? s.charCodeAt(c) : proceduralCode(r, c);
            count[k === 46 ? 45 : k]++;
          }
        }
        let best = 0;
        for (let k = 1; k < 128; k++) if (k !== 45 && count[k] > count[best]) best = k;
        if (best === 0) best = 45;
        out[sy * TILE + sx] = best;
        out[TILE * TILE + sy * TILE + sx] = Math.round((count[best] / n) * 255);
      }
    }
    return out;
  };
  const oneD = (n: number, lx: number, tx: number, f: (i: number) => number) => {
    const out = new Float64Array(2 * TILE).fill(Number.NaN);
    for (let s = 0; s < TILE; s++) {
      const from = (tx * TILE + s) * 2 ** lx;
      for (let i = from; i < Math.min(from + 2 ** lx, n); i++) {
        const v = f(i);
        if (!(out[2 * s] <= v)) out[2 * s] = v;
        if (!(out[2 * s + 1] >= v)) out[2 * s + 1] = v;
      }
    }
    return new Uint8Array(out.buffer);
  };
  const transport: TileTransport = async (path) => {
    calls.push(path);
    const [kind, ...rest] = path.split('/');
    const a = rest.map(Number);
    if (kind === 'tiles') return tile(a[0], a[1], a[2], a[3]);
    if (kind === 'labels') return json(Array.from({ length: a[1] - a[0] }, (_, i) => label(a[0] + i)));
    if (kind === 'rows' && a[0] < o.rowCount) return json({ label: label(a[0]), row: text(a[0]) });
    const stat = kind === 'stats' ? o.stats?.[rest[0]] : undefined;
    if (stat) return oneD(o.columnCount, a[1], a[2], stat);
    if (kind === 'values' && o.values) return oneD(o.length ?? 0, a[0], a[1], o.values);
    if (kind === 'counts' && o.counts) {
      const out = new Uint32Array(TILE);
      for (let s = 0; s < TILE; s++) {
        const from = (a[1] * TILE + s) * 2 ** a[0];
        for (let i = from; i < Math.min(from + 2 ** a[0], o.length ?? 0); i++) out[s] += o.counts(i);
      }
      return new Uint8Array(out.buffer);
    }
    throw new NotServed(`404 ${path}`);
  };
  return { transport, calls };
}
