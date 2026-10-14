/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tooltip follows the pointer without allocating (UX2): moves inside one hovered element change only the pointer,
// so each step is the tooltip's placement from offsets, and the keyboard cursor's anchor likewise. The
// viewer is loaded afresh for the case (`vi.resetModules`), as the surface's guards do (heap-browser.ts). The positive
// control: the moves allocating garbage each, sized to the run (`leaking`), warmed and measured alike, are caught.

import { expect, it, vi } from "vitest";
import { controlCaught, heapScaling, isFlat, leaking } from "../canvas/heap-browser.ts";

const N = 2_000;

/** The guard's measurement: 3,000 steps warm, then `N` and `4N`; and the next step's index. */
function measure(step: (i: number) => void) {
  for (let i = 0; i < 3_000; i++) step(i);
  let i = 0;
  const r = heapScaling(N, () => step(i++));
  return { r, next: i };
}

it("pointer moves inside one element place the tooltip in a flat heap", async () => {
  vi.resetModules();
  const { tooltipFixture } = await import("./tooltip-fixtures.ts");
  const t = await tooltipFixture();
  try {
    const { r, next } = measure(t.moveStep);
    let i = next;
    expect(isFlat(r), JSON.stringify(r)).toBe(true);
    // The steps placed it: it moved with the pointer, over one hovered element.
    const seen = new Set<string>();
    for (let k = 0; k < 8; k++) {
      t.moveStep(i++);
      seen.add(t.tip.style.getPropertyValue("--seqv-x"));
    }
    expect(seen.size).toBeGreaterThan(3);
    expect(t.hovers()).toBe(0);
  } finally {
    t.dispose();
  }
}, 120_000);

it("the keyboard cursor's anchor is placed in a flat heap", async () => {
  vi.resetModules();
  const { tooltipFixture } = await import("./tooltip-fixtures.ts");
  const t = await tooltipFixture();
  try {
    const { r } = measure(t.cursorStep());
    expect(isFlat(r), JSON.stringify(r)).toBe(true);
  } finally {
    t.dispose();
  }
}, 120_000);

it("the same moves allocating garbage sized to the run each are caught (positive control)", async () => {
  vi.resetModules();
  const { tooltipFixture } = await import("./tooltip-fixtures.ts");
  const t = await tooltipFixture();
  try {
    expect(controlCaught(() => measure(leaking(t.moveStep, N)).r)).toBe(true);
  } finally {
    t.dispose();
  }
}, 120_000);
