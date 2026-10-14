/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The self-test of heapGrowth: it reports a small allocation, and refuses a run that collects. And of heapScaling
// (RB0b): it catches a small allocation per step that one run under LIMIT does not, and passes a step that allocates
// nothing or allocates once per run.

import { describe, expect, it } from "vitest";
import {
  controlCaught,
  describeScaling,
  heapGrowth,
  heapScaling,
  isFlat,
  LIMIT,
  leaking,
  SCALED_LIMIT,
} from "./heap.ts";

const sink: { v: unknown } = { v: null };
const clean = () => {
  let s = 0;
  for (let i = 0; i < 1_000_000; i++) s += i & 7;
  sink.v = s;
};
/** `count` short-lived objects, none kept: 40 bytes each in Node (three header words and two fields, uncompressed). */
const garbage = (count: number) => () => {
  for (let i = 0; i < count; i++) sink.v = { a: i, b: i };
};

describe("heapGrowth", () => {
  it("reports no growth for a run that allocates nothing", () => {
    expect(heapGrowth(clean)).toBeLessThan(LIMIT);
  });

  it("reports a deliberate allocation that a single young generation holds", () => {
    expect(heapGrowth(garbage(200_000))).toBeGreaterThan(4 << 20); // 8 MB
  });

  it("fails a run that collects, though the heap ends small", () => {
    // 80 MB over the 64 MB young generation: one scavenge, after which what is left to see depends on where it fell
    // (1.7 million objects left 0.9 MB, 2026-10-11: its blind spot).
    expect(() => heapGrowth(garbage(2_000_000))).toThrow(/GC\(s\) ran during the measured run/);
  });
});

describe("heapScaling", () => {
  const n = 1000;
  const f32 = new Float32Array(256);

  it("passes a step that allocates nothing", () => {
    let s = 0;
    const r = heapScaling(n, (i) => {
      s = (s + (i & 7)) | 0;
      sink.v = s;
    });
    expect(isFlat(r), describeScaling(r, n)).toBe(true);
  });

  it("fails BK1's case, 208 bytes a step, that one run under LIMIT passes", () => {
    // Two `subarray` views a step, 104 bytes each (BK1's `PathBatch.append`, fixed there): 812 KB over 4,000 steps.
    const r = heapScaling(n, (i) => {
      sink.v = f32.subarray(i & 7);
      sink.v = f32.subarray(0, 8 + (i & 7));
    });
    expect(r.long, describeScaling(r, n)).toBeLessThan(LIMIT);
    expect(r.scaled).toBeGreaterThan(3 * n * 200);
    expect(isFlat(r)).toBe(false);
  });

  it("passes a step that allocates 512 KB once per run, not per step", () => {
    const r = heapScaling(n, (i) => {
      if (i === 0) sink.v = new Array(1 << 16).fill(i); // On the heap, as a typed array's store is not.
    });
    expect(r.short, describeScaling(r, n)).toBeGreaterThan(SCALED_LIMIT);
    expect(isFlat(r), describeScaling(r, n)).toBe(true);
  });

  it("refuses a run that collects", () => {
    // 40 KB a step, 160 MB over the longer run's 4,000 steps: past the 64 MB young generation.
    expect(() => heapScaling(n, () => garbage(1000)())).toThrow(/GC\(s\) ran during the measured run/);
  });
});

describe("the guard template's control", () => {
  const n = 500;
  let s = 0;
  const step = (i: number) => {
    s = (s + (i & 7)) | 0;
    sink.v = s;
  };

  it("catches the leaking step through a guard's own measurement, and passes the step alone", () => {
    expect(controlCaught(() => heapScaling(n, step))).toBe(false);
    expect(controlCaught(() => heapScaling(n, leaking(step)))).toBe(true);
  });

  it("counts a run refused for a collection as caught, as a single run's growth", () => {
    expect(controlCaught(() => heapScaling(n, () => garbage(1000)()))).toBe(true);
    expect(controlCaught(() => heapGrowth(garbage(200_000)))).toBe(true);
  });
});
