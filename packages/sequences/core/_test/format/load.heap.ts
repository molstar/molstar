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
// copied the rows). Untimed: `heapScaling` and `heapGrowth` run their own warm-ups, by count. The cases are
// load-cases.ts's, shared with the quick twin (load.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapGrowth, heapScaling, isFlat, leak, leaking } from '../model/heap.js';
import { delimited, STOCKHOLM_ROWS, STOCKHOLM_WIDTH, stockholm, TSV_ROWS, tileStep } from './load-cases.js';

describe('tile addressing', () => {
  it('allocates nothing', () => {
    const { step, sum } = tileStep();
    const n = 100_000;
    const r = heapScaling(n, step);
    expect(sum()).toBeGreaterThan(0);
    expect(isFlat(r), describeScaling(r, n)).toBe(true);
    // The template's control (CF-05): the same step with a leak is caught.
    expect(controlCaught(() => heapScaling(n, leaking(step)))).toBe(true);
  });
});

describe('a Stockholm tag table', () => {
  const [rows, width] = [STOCKHOLM_ROWS, STOCKHOLM_WIDTH];
  const file = stockholm();

  // The columns: values a double and a string pointer per position, 16 bytes; ranges (one per position here, each
  // digit differing from the next) two indices and two string pointers. The old reader's objects made it several times that.
  it.each([
    ['values', 20],
    ['ranges', 28],
  ] as const)('holds, for 600,000 positions as %s, little beyond the columns it returns', (geometry, limit) => {
    const table = file.table(geometry);
    let kept = 0;
    // The parse of the file is the model's; this is the table on top of it, so the parse's own garbage is subtracted.
    const parse = heapGrowth(() => {
      kept = file.parse() ? 1 : 0;
    });
    const total = heapGrowth(() => {
      kept = table() ? 1 : 0;
    });
    expect(kept).toBe(1);
    const perPosition = (total - parse) / (rows * width);
    expect(perPosition).toBeLessThan(limit);
    // The control (CF-05): the old reader's fault, an object per position, is caught by the same measurement and
    // bound, per position (LA-24).
    const objectEach = () => {
      table();
      for (let k = 0; k < rows * width; k++) leak(k);
    };
    expect(controlCaught(() => (heapGrowth(objectEach) - parse) / (rows * width), limit)).toBe(true);
  });
});

describe('a delimited file', () => {
  const rows = TSV_ROWS;
  const file = delimited();

  it('costs its cells: no array per row, no string built a character at a time', () => {
    let kept = 0;
    const decode = heapGrowth(() => {
      kept = file.decode();
    });
    const total = heapGrowth(() => {
      kept = file.read() ? 1 : 0;
    });
    expect(kept).toBe(1);
    // A cell's string (a short one is copied: a header and its characters) and its slot in the column.
    const perCell = (total - decode) / (rows * 4);
    expect(perCell).toBeLessThan(40);
    // The control (CF-05): an object per cell, garbage as the old reader's characters and row arrays were, is caught
    // by the same measurement and bound, per cell (LA-24).
    const objectEach = () => {
      file.read();
      for (let k = 0; k < rows * 4; k++) leak(k);
    };
    expect(controlCaught(() => (heapGrowth(objectEach) - decode) / (rows * 4), 40)).toBe(true);
  });
});
