/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A level built across frames (performance.md, Levels): the pyramid's build by a budget of cells a call gives the
// level its whole build gives; the rows keep the level drawn last while the view moves and the wanted one is not whole,
// and a settled frame builds the rest.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import { Pyramid } from '../pyramid.js';
import { rowsLevel } from '../rows-levels.js';

function base(rows: number, columns: number, seed: number): Uint8Array {
  const rng = new PCG(seed);
  const out = new Uint8Array(rows * columns);
  for (let k = 0; k < out.length; k++) out[k] = rng.int() % 5;
  return out;
}

function majority(rows: number, columns: number, seed = 1): Pyramid {
  const p = new Pyramid({ factor: 4, pairing: 'independent', reducer: 'majority' });
  p.setClasses(rows, columns, base(rows, columns, seed));
  return p;
}

describe('a pyramid level built across calls', () => {
  for (const [r, c] of [
    [2, 0],
    [1, 2],
    [4, 1],
  ] as const)
    it(`level (${r}, ${c}) by 5,000 cells a call equals its whole build`, () => {
      const whole = majority(3000, 200).level(r, c);
      const p = majority(3000, 200);
      let calls = 0;
      let l = p.levelWithin(r, c, 5000);
      while (!l) {
        calls++;
        expect(p.isBuilt(r, c)).toBe(false);
        l = p.levelWithin(r, c, 5000);
      }
      // a level of 600,000 base cells (or its finer levels') takes several calls
      expect(calls).toBeGreaterThan(3);
      expect([l.rows, l.columns]).toEqual([whole.rows, whole.columns]);
      const n = l.rows * l.columns;
      expect(Array.from(l.classes.subarray(0, n))).toEqual(Array.from(whole.classes.subarray(0, n)));
      expect(Array.from(l.shares.subarray(0, n))).toEqual(Array.from(whole.shares.subarray(0, n)));
    });

  it('a mean pyramid too; a new base builds again', () => {
    const values = new Float32Array(500 * 40).map((_, k) => (k * 7) % 13);
    const whole = new Pyramid({ factor: 2, pairing: 'square', reducer: 'mean' });
    whole.setValues(500, 40, values);
    const p = new Pyramid({ factor: 2, pairing: 'square', reducer: 'mean' });
    p.setValues(500, 40, values);
    let l = p.levelWithin(3, 3, 100);
    while (!l) l = p.levelWithin(3, 3, 100);
    expect(Array.from(l.values.subarray(0, l.rows * l.columns))).toEqual(
      Array.from(whole.level(3).values.subarray(0, l.rows * l.columns)),
    );
    p.setValues(500, 40, values);
    expect(p.isBuilt(3)).toBe(false);
  });
});

describe('the rows’ level while the view moves', () => {
  it('keeps the level drawn last until the wanted one is whole; a settled frame builds it', () => {
    const p = majority(3000, 200);
    // none drawn yet while the view moves: nothing until the level is whole (never a whole build in a moving frame)
    let first = rowsLevel(p, 1, 1, false, 5000);
    expect(first).toBeUndefined();
    while (!first) first = rowsLevel(p, 1, 1, false, 5000);
    expect(first).toMatchObject({ rowLevel: 1, columnLevel: 1 });
    // a zoom out: the old level drawn while the new one builds, a budget a frame
    let frames = 0;
    let l = rowsLevel(p, 3, 2, false, 5000);
    while (l!.rowLevel !== 3) {
      expect([l!.rowLevel, l!.columnLevel]).toEqual([1, 1]);
      frames++;
      l = rowsLevel(p, 3, 2, false, 5000);
    }
    expect(frames).toBeGreaterThan(1);
    // settled: built at once
    expect(rowsLevel(p, 0, 3, true, 5000)).toMatchObject({ rowLevel: 0, columnLevel: 3 });
  });
});
