/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The load paths' cases, for their heap guard and its quick twin (load.heap.ts, load.quick.ts): tile addressing (a step
// is a handful of `Tiles` calls, which the region cache makes per tile per frame), a Stockholm file read for its `GR:PP`
// tag table as values or as ranges, and a delimited file read whole. A read is one run: the quick twin measures the
// bytes a run allocates, the guard the heap it grows by.

import { TILE, Tiles } from '../../src/format/index.js';
import { KernelFormats } from '../../src/format/load/formats.js';

/** The tile-addressing step. The sum goes to a typed array's slot, never a captured `let`: past the small-integer
 * range (a few steps in), each store to the closure's context boxed a new heap number, 16 bytes, twice a step: the
 * guard's own 32 bytes a step. */
export function tileStep() {
  const sink = new Float64Array(1);
  const step = (i: number) => {
    const level = i & 7;
    sink[0] += Tiles.tileOf(i, level) + Tiles.start(i & 63, level) + Tiles.sampleOf(i, level);
    sink[0] += Tiles.count(i, level) + Tiles.levelFor(i + 1) + (i % TILE);
  };
  return { step, sum: () => sink[0] };
}

/** The reader's context: nothing reports. */
export const ctx = { report: () => {} } as never;

export const STOCKHOLM_ROWS = 2000;
export const STOCKHOLM_WIDTH = 300;

/** A Stockholm file of `rows` rows (the guard's 2,000) of 300 residues, with a `PP` tag line per row. */
export function stockholm(rows = STOCKHOLM_ROWS) {
  const width = STOCKHOLM_WIDTH;
  const names = Array.from({ length: rows }, (_, k) => `r${k}`);
  const letters = (k: number) => Array.from({ length: width }, (_, c) => 'ACDE'[(k + c) % 4]).join('');
  const digits = (k: number) => Array.from({ length: width }, (_, c) => ((k + c) % 10).toString()).join('');
  const data = new TextEncoder().encode(
    `# STOCKHOLM 1.0\n${names.map((n, k) => `${n} ${letters(k)}`).join('\n')}\n${names
      .map((n, k) => `#=GR ${n} PP ${digits(k)}`)
      .join('\n')}\n//\n`,
  );
  const format = KernelFormats.find((f) => f.name === 'stockholm');
  /** The file read without its tag table: the parse alone. */
  const parse = () => format?.read(data, { record: null, tag: null, geometry: null } as never, ctx);
  /** The tag table as `geometry`, on top of the parse. */
  const table = (geometry: 'values' | 'ranges') => () =>
    format?.read(data, { record: null, tag: 'GR:PP', geometry } as never, ctx);
  return { parse, table };
}

export const TSV_ROWS = 40_000;

/** A tab-separated file of 40,000 rows of four cells. */
export function delimited() {
  const lines = ['target\tstart\tstop\tlabel'];
  for (let i = 0; i < TSV_ROWS; i++) lines.push(`chr${i % 24}\t${i * 10}\t${i * 10 + 5}\t"site ${i % 7}"`);
  const data = new TextEncoder().encode(lines.join('\n'));
  const format = KernelFormats.find((f) => f.name === 'tsv');
  return {
    data,
    /** The read of the file. */
    read: () => format?.read(data, {} as never, ctx),
    /** The decoding of its bytes alone: what the read cannot avoid. */
    decode: () => new TextDecoder().decode(data).length,
  };
}
