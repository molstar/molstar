/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The React UI's pointer-rate and frame-rate paths allocate nothing (UX3): pointer moves inside one element place the
// tooltip and ask the track's box with no render, and a scroll of a few positions within a row moves the row controls
// by a custom property from cached strings. The viewer is loaded afresh per case (`vi.resetModules`), as the surface's
// guards do (heap-browser.ts). The positive control: the moves allocating garbage each, sized to the run (`leaking`),
// warmed and measured alike, are caught.

import { expect, it, vi } from "vitest";
import { controlCaught, heapScaling, isFlat, leaking } from "../canvas/heap-browser.js";

const N = 2_000;
// A scroll step is cheap and its string would be small (about 30 bytes): enough steps that one per step passes the limit.
const SCROLLS = 20_000;

/** The moves' measurement: 3,000 steps warm, `warmed` called, then `N` and `4N`; and the next step's index. */
function measure(step: (i: number) => void, warmed = () => {}) {
  for (let i = 0; i < 3_000; i++) step(i);
  warmed();
  let i = 0;
  const r = heapScaling(N, () => step(i++));
  return { r, next: i };
}

it("pointer moves inside one element, on the React UI, render nothing and allocate nothing", async () => {
  vi.resetModules();
  const { reactFixture } = await import("./heap-fixtures.js");
  const t = await reactFixture();
  try {
    let renders = 0;
    const { r, next } = measure(t.moveStep, () => {
      renders = t.renders();
    });
    let i = next;
    expect(isFlat(r), JSON.stringify(r)).toBe(true);
    expect(t.renders()).toBe(renders);
    const tip = t.host.querySelector(".seqv-tooltip") as HTMLElement;
    const seen = new Set<string>();
    for (let k = 0; k < 8; k++) {
      t.moveStep(i++);
      seen.add(tip.style.getPropertyValue("--seqv-x"));
    }
    expect(seen.size).toBeGreaterThan(3);
  } finally {
    await t.dispose();
  }
}, 120_000);

it("a scroll over a few positions moves the row controls in a flat heap", async () => {
  vi.resetModules();
  const { reactFixture } = await import("./heap-fixtures.js");
  const t = await reactFixture("msa");
  try {
    const inner = t.host.querySelector(".seqv-row-controls > div") as HTMLElement;
    expect(inner).not.toBeNull();
    for (let i = 0; i < 3_000; i++) t.scrollStep(i);
    const renders = t.renders();
    let i = 0;
    const r = heapScaling(SCROLLS, () => t.scrollStep(i++));
    expect(isFlat(r), JSON.stringify(r)).toBe(true);
    expect(t.renders()).toBe(renders);
    t.scrollStep(7);
    expect(inner.style.getPropertyValue("--seqv-scroll-y")).toBe(`${-(7 * 0.25)}px`);
  } finally {
    await t.dispose();
  }
}, 120_000);

it("a scroll far down a tall alignment (past 8,192 px) moves the row controls in a flat heap too", async () => {
  vi.resetModules();
  const { reactFixture } = await import("./heap-fixtures.js");
  const t = await reactFixture("msa");
  try {
    const inner = t.host.querySelector(".seqv-row-controls > div") as HTMLElement;
    // Across a page's edge: 20,479 px and on.
    const far = 20_479;
    for (let i = 0; i < 3_000; i++) t.scrollStep(i, far);
    let i = 0;
    const r = heapScaling(SCROLLS, () => t.scrollStep(i++, far));
    expect(isFlat(r), JSON.stringify(r)).toBe(true);
    t.scrollStep(7, far);
    const y = (name: string) => Number.parseFloat(inner.style.getPropertyValue(name));
    expect(y("--seqv-scroll-page") + y("--seqv-scroll-y")).toBe(-(far + 7 * 0.25));
  } finally {
    await t.dispose();
  }
}, 120_000);

it("the same moves allocating garbage sized to the run each are caught (positive control)", async () => {
  vi.resetModules();
  const { reactFixture } = await import("./heap-fixtures.js");
  const t = await reactFixture();
  try {
    expect(controlCaught(() => measure(leaking(t.moveStep, N)).r)).toBe(true);
  } finally {
    await t.dispose();
  }
}, 120_000);
