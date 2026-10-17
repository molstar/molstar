/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

import { afterAll, afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import type { HoverLine } from '../../../model/playback.js';
import { observableState } from '../../../reactive/observable-state.js';
import { classes } from '../_generated/hover-label.module.css.js';
import { HoverLabel } from '../hover-label.js';

// The hover label's frames (David, 2026-10-17): it appears at once (in the corner; beside the pointer in the frame that
// places it), and a hide waits one frame (a fast move between atoms, or out and back, does not blink it) and happens in
// that frame, fading out (`out`) with its last text, out of the accessibility tree.

beforeAll(() => {
  vi.useFakeTimers();
});
afterAll(() => {
  vi.useRealTimers();
});
afterEach(() => {
  document.body.innerHTML = '';
});
const frame = () => vi.advanceTimersByTime(17);

async function mount(placement: 'corner' | 'cursor') {
  const area = document.createElement('div');
  document.body.append(area);
  const lines = observableState<readonly HoverLine[]>(['ALA 12']);
  render(() => <HoverLabel model={{ lines: lines.asReadonly }} area={() => area} placement={placement} />, area);
  // The pointer's listeners attach once the label is in the page (a microtask).
  await Promise.resolve();
  const label = area.querySelector<HTMLElement>('[role=tooltip]')!;
  const move = () => area.dispatchEvent(new PointerEvent('pointermove', { clientX: 40, clientY: 30, bubbles: true }));
  const leave = () => area.dispatchEvent(new PointerEvent('pointerleave'));
  const shown = () => !label.classList.contains(classes.out) && !label.hasAttribute('aria-hidden');
  return { label, lines, move, leave, shown };
}

describe('hover label frames', () => {
  it('in the corner, shows at once when the pointer comes', async () => {
    const { move, shown } = await mount('corner');
    move();
    expect(shown()).toBe(true);
  });

  it('beside the pointer, shows in the frame that places it', async () => {
    const { move, shown } = await mount('cursor');
    move();
    frame();
    expect(shown()).toBe(true);
  });

  for (const placement of ['corner', 'cursor'] as const) {
    it(`${placement}: emptied, hides one frame later, in that frame, keeping its text to fade`, async () => {
      const { label, lines, move, shown } = await mount(placement);
      move();
      frame();
      expect(shown()).toBe(true);
      lines.set([]);
      expect(shown()).toBe(true);
      frame();
      expect(shown()).toBe(false);
      expect(label.getAttribute('aria-hidden')).toBe('true');
      expect(label.textContent).toBe('ALA 12');
    });

    it(`${placement}: the pointer gone, hides one frame later; back within the frame, stays`, async () => {
      const { move, leave, shown } = await mount(placement);
      move();
      frame();
      leave();
      move();
      frame();
      expect(shown()).toBe(true);
      leave();
      expect(shown()).toBe(true);
      frame();
      expect(shown()).toBe(false);
    });
  }
});
