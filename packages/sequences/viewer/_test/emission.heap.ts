/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for emission (seqspec:viewer#performance): a step a set, a bump or an emission, over runs of n and 4n
// (`heapScaling`, RB0b), each path warmed up first. The positive control: the state sets allocating an object a step
// (`leaking`), measured alike, are caught.

import { describe, expect, it } from "vitest";
import { event, live, state } from "../src/state.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.ts";

const N = 20_000;
const MASK = 0xfffff;

let sink = 0;

/** A state with two rxjs subscribers, and its step: one set. */
function setting() {
  const s = state(0);
  s.subscribe((v) => {
    sink = (sink + v) & MASK;
  });
  s.subscribe((v) => {
    sink = (sink ^ v) & MASK;
  });
  return (k: number) => s.set(k & MASK, "pointer");
}

describe("emission allocates nothing", () => {
  it("state sets, with rxjs subscribers", () => {
    const r = heapScaling(N, setting());
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("state sets queued during delivery", () => {
    const s = state(0);
    // An odd value sets the next one from inside the delivery: each outer set queues one.
    s.subscribe((v) => {
      if (v & 1) s.set((v + 1) & MASK, "behavior");
    });
    s.subscribe((v) => {
      sink = (sink + v) & MASK;
    });
    const r = heapScaling(N, (k) => s.set((2 * k + 1) & MASK, "pointer"));
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("live version bumps", () => {
    const l = live({ x: 0 });
    l.version.subscribe((v) => {
      sink = (sink + v + l.current.x) & MASK;
    });
    const r = heapScaling(N, (k) => {
      l.current.x = k & MASK;
      l.bump();
    });
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("events with a reused record, by priority", () => {
    const e = event<{ x: number; owner: string }>();
    const rec = { x: 0, owner: "" };
    e.subscribe(
      (r) => {
        sink = (sink + r.x) & MASK;
      },
      { priority: 100 },
    );
    e.subscribe(
      (r) => {
        if (r.x === -1) r.owner = "never";
      },
      { priority: 50 },
    );
    e.subscribe((r) => {
      sink = (sink ^ r.x) & MASK;
    });
    const r = heapScaling(N, (k) => {
      rec.x = k & MASK;
      e.emit(rec);
    });
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(sink).toBeGreaterThanOrEqual(0);
  });

  // The control last: cases in one file share V8's feedback.
  it("the same sets allocating an object each are caught (positive control)", () => {
    expect(controlCaught(() => heapScaling(N, leaking(setting())))).toBe(true);
  });
});
