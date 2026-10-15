/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of keyboard.heap.ts (`test:perf`): IN4's key repeat (keyboard-cases.ts, the guard's cases), each with
// its control. A press that leaves the cursor reads nothing; one that moves it, the one record the state publishes (the
// guard's per-press bound, 112 bytes).

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';

/** The cursor moves' counts: they allocate by design, never read quiet, and run to this cap. */
const BY_DESIGN = { maxWarmSteps: 6_000 };
import { along, family, millionRows, rects, staying, toBottomEnd, upDown, zooming } from './keyboard-cases.js';

describe("IN4's key repeat, quick", () => {
  it('presses that leave the cursor: the end of a row, the bottom, the menu key', async ({ task }) => {
    const t = await family();
    toBottomEnd(t);
    const at = t.ctx.ui.keyCursor.value;
    quickCase(task, 'keyboard: staying', staying(t));
    expect(t.ctx.ui.keyCursor.value).toBe(at);
  });

  it('cursorRect on data, a row and a header', async ({ task }) => {
    quickCase(task, 'keyboard: cursorRect', rects(await family()));
  });

  it('zoom keys at the cursor', async ({ task }) => {
    const t = await family();
    const step = zooming(t);
    const at = t.ctx.ui.keyCursor.value;
    quickCase(task, 'keyboard: zoom keys', step);
    expect(t.ctx.ui.keyCursor.value).toBe(at);
  });

  it('Left and Right along a row', async ({ task }) => {
    // A cursor move allocates by design (80 bytes a press), so never quiet: a short cap (alloc.ts).
    quickCase(task, 'keyboard: along a row', along(await family()), BY_DESIGN);
  });

  it('Up and Down, Home and End', async ({ task }) => {
    quickCase(task, 'keyboard: up and down', upDown(await family()), BY_DESIGN);
  });

  it('Down and Up over a million-row run', async ({ task }) => {
    quickCase(task, 'keyboard: million rows', await millionRows(), BY_DESIGN);
  });
});
