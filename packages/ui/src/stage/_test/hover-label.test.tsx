/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { DomHost, render } from '../../dom/index.js';
import { classes } from '../_generated/hover-label.module.css.js';
import { HoverLabel } from '../hover-label.js';
import type { HoverLine } from '../stage-model.js';

// The hover label (stage §3.1): in the corner by default, beside the pointer when asked; text only; it appears at once
// (in the corner; beside the pointer in the frame that places it), a hide waits one frame and happens in it, fading
// with its last text, out of the accessibility tree; a held button hides it at once.

const host = DomHost.requestAnimationFrame;
let queue: ((t: number) => void)[] = [];
let now = 0;
const frame = () => {
  now += 16;
  const run = queue;
  queue = [];
  for (const f of run) f(now);
};

beforeEach(() => {
  DomHost.requestAnimationFrame = (f) => {
    queue.push(f);
    return 0;
  };
});
afterEach(() => {
  document.body.innerHTML = '';
  for (let i = 0; i < 4 && queue.length; i++) frame();
  DomHost.requestAnimationFrame = host;
  queue = [];
});

async function mount(placement: 'corner' | 'cursor', initial: readonly HoverLine[] = ['ALA 12']) {
  const area = document.createElement('div');
  document.body.append(area);
  const lines = new State<readonly HoverLine[]>(initial);
  const dispose = render(() => <HoverLabel model={{ lines }} area={() => area} placement={placement} />, area);
  // The pointer's listeners attach once the label is in the page (a microtask).
  await Promise.resolve();
  const label = area.querySelector<HTMLElement>('[role=tooltip]')!;
  const move = (buttons = 0) =>
    area.dispatchEvent(new PointerEvent('pointermove', { clientX: 40, clientY: 30, buttons, bubbles: true }));
  const leave = () => area.dispatchEvent(new PointerEvent('pointerleave'));
  const shown = () => !label.classList.contains(classes.out) && !label.hasAttribute('aria-hidden');
  return { label, lines, move, leave, shown, dispose };
}

describe('hover label placement', () => {
  it('is in the corner by default and never moves there', async () => {
    const { label, move, shown } = await mount('corner');
    expect(label.classList.contains(classes.corner)).toBe(true);
    move();
    expect(shown()).toBe(true);
    frame();
    expect(label.style.getPropertyValue('--msui-hover-x')).toBe('');
  });

  it('beside the pointer, shows in the frame that places it', async () => {
    const { label, move, shown } = await mount('cursor');
    expect(label.classList.contains(classes.cursor)).toBe(true);
    move();
    expect(shown()).toBe(false);
    frame();
    expect(shown()).toBe(true);
    // No layout in jsdom: the offset parent is none, so nothing is placed but it shows.
    expect(label.getAttribute('aria-hidden')).toBeNull();
  });
});

describe('hover label content', () => {
  it('renders runs as elements, never as markup', async () => {
    const { label } = await mount('corner', [
      { role: 'title', runs: [{ text: 'PHE', bold: true }, { text: ' <img src=x onerror=alert(1)>' }] },
      {
        role: 'text',
        runs: [
          { text: 'H', code: true },
          { text: '2', sub: true },
        ],
      },
    ]);
    expect(label.querySelector('strong')?.textContent).toBe('PHE');
    expect(label.querySelector('code')?.textContent).toBe('H');
    expect(label.querySelector('sub')?.textContent).toBe('2');
    expect(label.querySelector('img')).toBeNull();
    expect(label.textContent).toContain('<img src=x onerror=alert(1)>');
  });

  it('keeps the last label a frame when it empties, and swaps in the next without hiding', async () => {
    const { label, lines } = await mount('corner');
    lines.set([]);
    expect(label.textContent).toBe('ALA 12');
    lines.set(['GLY 13']);
    expect(label.textContent).toBe('GLY 13');
    frame();
    expect(label.textContent).toBe('GLY 13');
  });
});

describe('hover label frames', () => {
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

    it(`${placement}: hidden at once while a button is held`, async () => {
      const { move, shown } = await mount(placement);
      move();
      frame();
      expect(shown()).toBe(true);
      move(1);
      expect(shown()).toBe(false);
    });
  }

  it('asks for no frame at rest', async () => {
    const { move } = await mount('cursor');
    move();
    frame();
    frame();
    expect(queue.length).toBe(0);
  });

  it('its listeners go with its scope', async () => {
    const { move, dispose, label, shown } = await mount('cursor');
    dispose();
    move();
    expect(queue.length).toBe(0);
    expect(shown()).toBe(false);
    expect(label.isConnected).toBe(false);
  });
});
