/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick allocation tier (`test:perf`) of the glass wordmark (glass.ts): a step and a draw of the pieces moving
// allocate nothing (everything a frame touches is preallocated, the stroke and the atlas looked up only when the colors
// or the scale change).

import { quickCase } from '@molstar/test-harness/quick';
import { expect, it } from 'vitest';
import { Glass, type GlassColors } from '../glass.js';
import { fakeContext, stubPath2D } from './fake-canvas.js';

stubPath2D();
const colors: GlassColors = { logo: ['#111', '#222', '#333', '#444'], edge: '#999' };

const cases: Record<string, () => (i: number) => void> = {
  /** In pieces under a moving pointer: a step a frame. */
  step() {
    const g = new Glass(23, 0.35);
    g.shatter(20, 10);
    return (i) => {
      g.setPointer((i % 50) + 2, 10);
      g.step(1 / 60);
      if (g.phase === 'whole') g.shatter(20, 10);
    };
  },
  /** Whole, the pointer moving over it: a step a frame (the hover's lean only). */
  lean() {
    const g = new Glass(23, 0.35);
    return (i) => {
      g.setPointer((i % 50) + 2, 10);
      g.step(1 / 60);
    };
  },
  /** Apart (held there: no gathering), a step a frame. */
  apart() {
    const g = new Glass(23, 0.35);
    g.shatter(20, 10);
    return (i) => {
      if ((i & 31) === 0) g.shatter(20, 10);
      g.step(1 / 600);
    };
  },
  /** In pieces: a draw a frame (nothing moved). */
  draw() {
    const g = new Glass(23, 0.35);
    const { ctx } = fakeContext();
    g.shatter(20, 10);
    g.step(1 / 60);
    return () => g.draw(ctx, 0, 0, 1, colors);
  },
  /** In pieces under a moving pointer: a step and a draw a frame. */
  stepAndDraw() {
    const g = new Glass(23, 0.35);
    const { ctx } = fakeContext();
    g.shatter(20, 10);
    return (i) => {
      g.setPointer((i % 50) + 2, 10);
      g.step(1 / 60);
      if (g.phase === 'whole') g.shatter(20, 10);
      g.draw(ctx, 0, 0, 1, colors);
    };
  },
};

for (const [name, make] of Object.entries(cases)) {
  it(name, ({ task }) => {
    quickCase(task, `glass: ${name}`, make());
    expect(true).toBe(true);
  });
}
