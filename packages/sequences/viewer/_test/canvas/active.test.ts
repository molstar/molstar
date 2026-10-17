/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A paused surface (`setActive(false)`; seqspec:viewer#the-render-backend), in Chromium: an embedding host hides it, and
// it asks for no frame while pans, hovers and resizes reach it; the flags wait. Resumed, one frame draws all that waited,
// at the size it was given while paused. The positive control: the same steps on an active surface ask for frames.

import { afterEach, expect, it } from 'vitest';
import { canvasesOf, pixelsOf } from './open.js';
import { surfaceFrames } from './surface-fixtures.js';

let t: Awaited<ReturnType<typeof surfaceFrames>> | null = null;
afterEach(() => {
  t?.dispose();
  t = null;
});

/** The scheduler's frame requests over `n` steps, and the surface's frames drawn. */
function over(n: number) {
  const s = t!;
  const asked = s.scheduler.requested;
  const frames = s.surface.stats.frames;
  for (let i = 0; i < n; i++) s.step(i);
  return { asked: s.scheduler.requested - asked, drawn: s.surface.stats.frames - frames };
}

it('an active surface asks for a frame per pan step and draws it (positive control)', async () => {
  t = await surfaceFrames('both');
  const r = over(20);
  expect(r.asked).toBeGreaterThanOrEqual(20);
  expect(r.drawn).toBe(20);
});

it('paused, a surface asks for no frame and draws none; resumed, one frame draws what waited', async () => {
  t = await surfaceFrames('both');
  const s = t;
  s.surface.setActive(false);
  const before = pixelsOf(canvasesOf({ surface: s.surface } as never)[0]);
  // The steps' events still run (the viewport pans, the hover moves): only the surface's frames are gone.
  // Odd: the pans (±9 px in turn) leave the view moved.
  const paused = over(51);
  expect(paused.drawn).toBe(0);
  // The pan's own frames are the surface's: none asked for. (The dispatcher's re-hit after a move may ask for one.)
  const panOnly = await surfaceFrames('pan');
  try {
    panOnly.surface.setActive(false);
    const asked = panOnly.scheduler.requested;
    const frames = panOnly.surface.stats.frames;
    for (let i = 0; i < 50; i++) panOnly.step(i);
    expect(panOnly.scheduler.requested - asked).toBe(0);
    expect(panOnly.surface.stats.frames).toBe(frames);
  } finally {
    panOnly.dispose();
  }
  // Resumed: one frame asked for, drawing the panned view.
  const asked = s.scheduler.requested;
  const frames = s.surface.stats.frames;
  s.surface.setActive(true);
  expect(s.scheduler.requested - asked).toBe(1);
  s.scheduler.flush();
  expect(s.surface.stats.frames).toBe(frames + 1);
  const after = pixelsOf(canvasesOf({ surface: s.surface } as never)[0]);
  expect(after.length).toBe(before.length);
  expect(after.some((v, i) => v !== before[i])).toBe(true);
});

it('a resize while paused is recorded, not drawn; the resume draws at the new size', async () => {
  t = await surfaceFrames('pan');
  const s = t;
  s.surface.setActive(false);
  const frames = s.surface.stats.frames;
  const host = s.surface.element.parentElement as HTMLElement;
  host.style.width = '500px';
  host.style.height = '300px';
  // The ResizeObserver reports in the next rendering frame.
  for (let k = 0; k < 2; k++) await new Promise((r) => requestAnimationFrame(r));
  expect(s.surface.stats.frames).toBe(frames);
  s.surface.setActive(true);
  s.scheduler.flush();
  expect(s.surface.stats.frames).toBe(frames + 1);
  const [main] = canvasesOf({ surface: s.surface } as never);
  // The fixture's DPR is 2.
  expect([main.width, main.height]).toEqual([1000, 600]);
});
