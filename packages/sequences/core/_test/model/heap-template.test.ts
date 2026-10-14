/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The guard template's control logic, measured by stand-ins (the heap tier measures for real: heap-guard.heap.ts).
import { describe, expect, it } from 'vitest';
import { controlCaught, LIMIT, leaking, SCALED_LIMIT } from './heap.js';

describe('controlCaught', () => {
  it('is caught when not flat, a growth past the limit, or a run refused for a collection', () => {
    expect(controlCaught(() => ({ short: 0, long: 0, scaled: 0 }))).toBe(false);
    expect(controlCaught(() => ({ short: 0, long: SCALED_LIMIT, scaled: SCALED_LIMIT }))).toBe(true);
    expect(controlCaught(() => LIMIT - 1)).toBe(false);
    expect(controlCaught(() => LIMIT)).toBe(true);
    expect(
      controlCaught(() => {
        throw new Error('1 GC(s) ran during the measured run, twice: it hides allocation. Make the run shorter.');
      }),
    ).toBe(true);
  });

  it("rethrows the guard's own failure", () => {
    expect(() =>
      controlCaught(() => {
        throw new Error('no section');
      }),
    ).toThrow('no section');
  });

  it("leaking runs the guard's step, then the leak", () => {
    const seen: number[] = [];
    const step = leaking((i) => seen.push(i));
    step(3);
    step(4);
    expect(seen).toEqual([3, 4]);
  });
});
