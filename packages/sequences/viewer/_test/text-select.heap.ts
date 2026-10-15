/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A select drag in a wrapped section allocates nothing per pointer event or per frame after warm-up
// (seqspec:interaction#selection; seqspec:viewer#performance): on the text preset, through the real dispatcher and
// behaviors, a press on the letters, then moves across four lines, past their ends and over their labels bands, the
// overlay drawing the pending range each step, and over the features' underlines by lane (each hit an item); and the
// scroll past the rows' edges, a frame a step, the range following.
// The positive control: the same moves allocating an object each (`leaking`), warmed and measured alike, are caught.
// Each case loads the viewer afresh (`vi.resetModules`) and is warmed and settled by count (`warmFrame`); the cases are
// text-select-cases.ts's, shared with the quick twin (text-select.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking, warmFrame } from './heap.js';
import { textSelect } from './text-select-cases.js';

/** A deep guard's steps: 500 catch about 43 bytes a step, these about 22. */
const STEPS = 1_000;

/**
 * The guard's measurement of a drag's step: warmed and settled by count, counted anew, then `STEPS` and `4 * STEPS`,
 * `step` getting a running index throughout (the scroll's frame clock runs on it: restarted at 0 it would run back).
 */
const measure = (d: Awaited<ReturnType<typeof textSelect>>, step: (i: number) => void) => {
  let k = 0;
  const run = () => step(k++);
  warmFrame(run);
  d.resetCounts();
  return heapScaling(STEPS, run);
};

async function run(mode: 'lines' | 'scroll') {
  const d = await textSelect(mode);
  const r = measure(d, d.step);
  return { r, ...d.state() };
}

describe('a select drag in text allocates nothing after warm-up', () => {
  it('moves across four lines, past their ends and over their labels bands, the pending range drawn each step', async () => {
    const { r, owner, ends, active } = await run('lines');
    expect([owner, active]).toEqual(['select', true]);
    // The range followed the pointer: its end changed on most steps.
    expect(ends).toBeGreaterThan(STEPS);
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 300_000);

  it("the scroll past the rows' edges, a frame a step, the range following", async () => {
    const { r, owner, ends, scrolls } = await run('scroll');
    // The rows scrolled on every step but the turns, and the range's end followed them across lines.
    expect(owner).toBe('select');
    expect(scrolls).toBeGreaterThan(9 * STEPS);
    expect(ends).toBeGreaterThan(STEPS);
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 300_000);

  it('the same moves allocating an object each are caught (positive control)', async () => {
    const d = await textSelect('lines');
    expect(controlCaught(() => measure(d, leaking(d.step)))).toBe(true);
    expect(d.state().owner).toBe('select');
  }, 300_000);
});
