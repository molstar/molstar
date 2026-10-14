/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `metricsMeasurer` allocates nothing per call but the number it returns (layout prepares labels through it at edge rate, a
// long list at a time): `metrics` returns a record made once per font string; `measure` over strings of ASCII and of other
// code points makes no string, no array and no record, and V8 boxes the fractional width it returns to a caller that does
// not inline it (16 bytes), which is all the guard allows. The positive control: `metrics` calls allocating an object
// each (`leaking`), through the same measurement, are caught.

import { expect, it } from "vitest";
import { DefaultFontMetrics, type EmBox, metricsMeasurer } from "../../src/svg/index.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking, SCALED_LIMIT } from "../heap.ts";

const N = 100_000;
const m = metricsMeasurer(DefaultFontMetrics);
const fonts = ["400 12px ui-sans-serif, system-ui, sans-serif", "600 11px ui-monospace, Menlo, monospace"];
const texts = ["chr12:1,234,567-1,299,999", "Zinc ligands (UniProt) é", "ABCDEFGHIKLMNPQRSTVWY"];

it("metrics() allocates nothing once a font string is known", () => {
  const sum = new Float64Array(1);
  const step = (i: number) => {
    sum[0] += (m.metrics(fonts[i & 1]) as EmBox).ascent;
  };
  const r = heapScaling(N, step);
  expect(sum[0]).toBeGreaterThan(0);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
  expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
});

it("measure() allocates only the width it returns", () => {
  const sum = new Float64Array(1); // a captured `let` double would be boxed on every step
  const r = heapScaling(N, (i) => {
    sum[0] += m.measure(fonts[i & 1], texts[i % 3]);
  });
  expect(sum[0]).toBeGreaterThan(0);
  expect(r.scaled, describeScaling(r, N)).toBeLessThanOrEqual(3 * N * 16 + SCALED_LIMIT);
});
