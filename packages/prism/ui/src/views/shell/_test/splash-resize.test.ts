/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

import { describe, expect, it } from 'vitest';
import { PrismField } from '../prism-field.js';
import { QuantumField } from '../quantum-field.js';

// A resize wipes the canvas (its size set), and the resize observer that saw it runs past the frame's tick: drawn only
// in the next frame, the splash would show one blank frame each step of a window's resize (a blink). `redraw` draws it
// at once (prism-splash.tsx, `place`).

/** A 2D context that takes every call and counts the clears (a draw starts with one). */
function stubContext() {
  const counts = { clears: 0 };
  const any: object = new Proxy(() => any, {
    get: (_t, key) => (key === 'clearRect' ? () => counts.clears++ : key === Symbol.toPrimitive ? () => 0 : any),
    apply: () => any,
    set: () => true,
  });
  return { ctx: any, counts };
}

function canvasOf(width: number, height: number) {
  const canvas = document.createElement('canvas');
  const { ctx, counts } = stubContext();
  (canvas as unknown as { getContext: () => unknown }).getContext = () => ctx;
  const size = { width, height };
  Object.defineProperty(canvas, 'clientWidth', { get: () => size.width });
  Object.defineProperty(canvas, 'clientHeight', { get: () => size.height });
  return { canvas, counts, size };
}

describe('the splash on a resize', () => {
  for (const [name, Field] of [
    ['light', PrismField],
    ['quantum', QuantumField],
  ] as const) {
    it(`${name}: draws at once after the canvas is resized, not a frame later`, () => {
      const { canvas, counts, size } = canvasOf(400, 300);
      const field = new Field(canvas, {}, {});
      field.resize();
      field.start();
      size.width = 500;
      field.resize();
      const before = counts.clears;
      field.redraw();
      expect(counts.clears).toBe(before + 1);
      field.dispose();
    });
  }

  it('draws nothing before it starts (nothing was drawn to keep)', () => {
    const { canvas, counts } = canvasOf(400, 300);
    const field = new PrismField(canvas, {}, {});
    field.resize();
    field.redraw();
    expect(counts.clears).toBe(0);
    field.dispose();
  });
});
