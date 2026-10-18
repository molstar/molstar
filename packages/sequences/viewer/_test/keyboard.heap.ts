/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for IN4's key repeat (seqspec:viewer#behaviors; the cursor is a published value, `ui.keyCursor`): on
// `heapScaling` (a step is one key press; labels off, since a label is a string per hovered element by design, bounded
// by IN2's churn guard). A press that leaves the cursor where it is allocates nothing: the edges of a row and of the
// rows, the zoom keys, the menu key (asked at the cursor, with no `interaction.menu` subscriber:
// seqspec:interaction#streams). A press that moves it makes the one small record the state publishes and nothing else:
// a bound per press, whatever the number of positions, rows and blocks (a million-row run included). The cases are
// keyboard-cases.ts', shared with the quick twin (keyboard.quick.ts). The positive controls, one per criterion: the
// presses that leave the cursor allocating an object each (`leaking`), primed and measured alike, are caught as not
// flat, and the presses that move it, allocating likewise, by the per-press bound.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { along, family, millionRows, rects, staying, toBottomEnd, upDown, zooming } from './keyboard-cases.js';

const N = 20_000;
/** One published cursor (88 bytes; 93 measured, 2026-10-11) is all a press may make. */
const PER_PRESS = 112;

/** At least 300,000 presses, then until three chunks in a row allocate under 2 KB: V8's lower tiers box doubles. */
function prime(step: (i: number) => void) {
  const gc = (globalThis as { gc?: () => void }).gc;
  let k = 0;
  for (let chunk = 0, quiet = 0; chunk < 100 && (quiet < 3 || chunk < 15); chunk++) {
    gc?.();
    const before = process.memoryUsage().heapUsed;
    for (let j = 0; j < 20_000; j++) step(k++);
    quiet = process.memoryUsage().heapUsed - before < 2048 ? quiet + 1 : 0;
  }
}
const guard = (step: (i: number) => void) => {
  prime(step);
  const r = heapScaling(N, step);
  console.info(`key repeat: ${(r.scaled / (3 * N)).toFixed(1)} bytes a press`);
  return r;
};

describe("IN4's key repeat", { timeout: 300_000 }, () => {
  it('a press that leaves the cursor where it is allocates nothing: the end of a row, the bottom, the top, the menu key', async () => {
    const t = await family();
    toBottomEnd(t);
    const at = t.ctx.ui.keyCursor.value;
    const r = guard(staying(t));
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(t.ctx.ui.keyCursor.value).toBe(at);
  });

  it('cursorRect, read per frame by the ring and the tooltip, allocates nothing: on data, a row and a header', async () => {
    const r = guard(rects(await family()));
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it('zoom keys, anchored at the cursor, allocate nothing; the cursor is not republished', async () => {
    const t = await family();
    const step = zooming(t);
    const at = t.ctx.ui.keyCursor.value;
    const r = guard(step);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(t.ctx.ui.keyCursor.value).toBe(at);
  });

  it('Left and Right along a row, across the header and back: a bounded record per press, flat in the count', async () => {
    const r = guard(along(await family()));
    expect(r.scaled / (3 * N), describeScaling(r, N)).toBeLessThan(PER_PRESS);
  });

  it('Up and Down along the tracks and rows, and Home and End: the same bound', async () => {
    const r = guard(upDown(await family()));
    expect(r.scaled / (3 * N), describeScaling(r, N)).toBeLessThan(PER_PRESS);
  });

  it('the rows of a million-row run: Down and Up keep the bound, however far the cursor is from the first row', async () => {
    const r = guard(await millionRows());
    expect(r.scaled / (3 * N), describeScaling(r, N)).toBeLessThan(PER_PRESS);
  });

  it('the same presses allocating an object each are caught (positive control)', async () => {
    const t = await family();
    toBottomEnd(t);
    expect(controlCaught(() => guard(leaking(staying(t))))).toBe(true);
  });

  it('presses that move the cursor, allocating an object each, are caught by the per-press bound (positive control)', async () => {
    const t = await family();
    expect(controlCaught(() => guard(leaking(along(t))).scaled / (3 * N), PER_PRESS)).toBe(true);
  });
});
