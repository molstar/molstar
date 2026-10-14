/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the load paths (VW3). Addressing a tile allocates nothing (`Tiles`, called per tile per frame by the
// region cache) on RB0b's stepped `heapScaling`; and a one-shot build, a Stockholm tag table, keeps what it returns and
// little else (the old reader built an object per position: for this table its garbage passes the 64 MB young
// generation, so a GC runs in the measured run and `heapGrowth` refuses it); and a delimited file keeps its cells and
// makes nothing else (the reader before the review built each cell a character at a time and an array per row, then
// copied the rows).

import { describe, expect, it } from 'vitest';
import { TILE, Tiles } from '../../src/format/index.js';
import { KernelFormats } from '../../src/format/load/formats.js';
import { controlCaught, describeScaling, heapGrowth, heapScaling, isFlat, leak, leaking } from '../model/heap.js';

describe('tile addressing', () => {
  it('allocates nothing', () => {
    // The sum goes to a typed array's slot, never a captured `let`: past the small-integer range (a few steps in), each
    // store to the closure's context boxed a new heap number, 16 bytes, twice a step: the guard's own 32 bytes a step.
    const sink = new Float64Array(1);
    const n = 100_000;
    const step = (i: number) => {
      const level = i & 7;
      sink[0] += Tiles.tileOf(i, level) + Tiles.start(i & 63, level) + Tiles.sampleOf(i, level);
      sink[0] += Tiles.count(i, level) + Tiles.levelFor(i + 1) + (i % TILE);
    };
    const r = heapScaling(n, step);
    expect(sink[0]).toBeGreaterThan(0);
    expect(isFlat(r), describeScaling(r, n)).toBe(true);
    // The template's control (CF-05): the same step with a leak is caught.
    expect(controlCaught(() => heapScaling(n, leaking(step)))).toBe(true);
  });
});

describe('a Stockholm tag table', () => {
  const [rows, width] = [2000, 300];
  const names = Array.from({ length: rows }, (_, k) => `r${k}`);
  const letters = (k: number) => Array.from({ length: width }, (_, c) => 'ACDE'[(k + c) % 4]).join('');
  const digits = (k: number) => Array.from({ length: width }, (_, c) => ((k + c) % 10).toString()).join('');
  const data = new TextEncoder().encode(
    `# STOCKHOLM 1.0\n${names.map((n, k) => `${n} ${letters(k)}`).join('\n')}\n${names
      .map((n, k) => `#=GR ${n} PP ${digits(k)}`)
      .join('\n')}\n//\n`,
  );
  const format = KernelFormats.find((f) => f.name === 'stockholm');

  // The columns: values a double and a string pointer per position, 16 bytes; ranges (one per position here, each
  // digit differing from the next) two indices and two string pointers. The old reader's objects made it several times that.
  it.each([
    ['values', 20],
    ['ranges', 28],
  ] as const)('holds, for 600,000 positions as %s, little beyond the columns it returns', (geometry, limit) => {
    const ctx = { report: () => {} } as never;
    const request = { record: null, tag: 'GR:PP', geometry } as never;
    let kept = 0;
    // The parse of the file is the model's; this is the table on top of it, so the parse's own garbage is subtracted.
    const parse = heapGrowth(() => {
      kept = format?.read(data, { ...(request as object), tag: null, geometry: null } as never, ctx) ? 1 : 0;
    });
    const total = heapGrowth(() => {
      kept = format?.read(data, request, ctx) ? 1 : 0;
    });
    expect(kept).toBe(1);
    const perPosition = (total - parse) / (rows * width);
    expect(perPosition).toBeLessThan(limit);
    // The control (CF-05): the old reader's fault, an object per position, is caught by the same measurement and
    // bound, per position (LA-24).
    const objectEach = () => {
      format?.read(data, request, ctx);
      for (let k = 0; k < rows * width; k++) leak(k);
    };
    expect(controlCaught(() => (heapGrowth(objectEach) - parse) / (rows * width), limit)).toBe(true);
  });
});

describe('a delimited file', () => {
  const rows = 40_000;
  const lines = ['target\tstart\tstop\tlabel'];
  for (let i = 0; i < rows; i++) lines.push(`chr${i % 24}\t${i * 10}\t${i * 10 + 5}\t"site ${i % 7}"`);
  const data = new TextEncoder().encode(lines.join('\n'));
  const format = KernelFormats.find((f) => f.name === 'tsv');

  it('costs its cells: no array per row, no string built a character at a time', () => {
    const ctx = { report: () => {} } as never;
    let kept = 0;
    const decode = heapGrowth(() => {
      kept = new TextDecoder().decode(data).length;
    });
    const total = heapGrowth(() => {
      kept = format?.read(data, {} as never, ctx) ? 1 : 0;
    });
    expect(kept).toBe(1);
    // A cell's string (a short one is copied: a header and its characters) and its slot in the column.
    const perCell = (total - decode) / (rows * 4);
    expect(perCell).toBeLessThan(40);
    // The control (CF-05): an object per cell, garbage as the old reader's characters and row arrays were, is caught
    // by the same measurement and bound, per cell (LA-24).
    const objectEach = () => {
      format?.read(data, {} as never, ctx);
      for (let k = 0; k < rows * 4; k++) leak(k);
    };
    expect(controlCaught(() => (heapGrowth(objectEach) - decode) / (rows * 4), 40)).toBe(true);
  });
});
