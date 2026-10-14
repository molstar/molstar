/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The self-test of the browser guard's measurement: a flat step passes, a step that allocates 1 KB does not, and a run
// that collects is refused.

import { expect, it } from "vitest";
import { controlCaught, heapScaling, isFlat, leaking } from "./heap-browser.js";

const sink: { v: unknown } = { v: null };

it("a step that allocates nothing is flat", () => {
  let s = 0;
  expect(
    isFlat(
      heapScaling(20_000, (i) => {
        s += i & 7;
        sink.v = s;
      }),
    ),
  ).toBe(true);
});

it("a step that allocates 1 KB scales with the run length", () => {
  const r = heapScaling(8_000, (i) => {
    sink.v = new Float64Array(128).fill(i);
  });
  expect(r.scaled).toBeGreaterThan(1 << 20);
  expect(isFlat(r)).toBe(false);
});

it("refuses a run that collects, though both lengths could end under the limit", () => {
  // About 100 MB a run, over the 64 MB young generation; near 64 MB both runs' growth wraps under 1 MB.
  expect(() =>
    heapScaling(100_000, (i) => {
      sink.v = new Array(256).fill(i);
    }),
  ).toThrow(/GC ran during a measured run/);
});

it("the guard template's control: a flat step leaking is caught, the step alone is not", () => {
  let s = 0;
  const step = (i: number) => {
    s += i & 7;
    sink.v = s;
  };
  expect(controlCaught(() => heapScaling(20_000, step))).toBe(false);
  expect(controlCaught(() => heapScaling(20_000, leaking(step)))).toBe(true);
});
