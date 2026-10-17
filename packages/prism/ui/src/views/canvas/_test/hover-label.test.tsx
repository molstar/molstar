/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

import { afterEach, describe, expect, it } from 'vitest';
import { render } from '../../../dom/render.js';
import type { HoverLine } from '../../../model/playback.js';
import { observableState } from '../../../reactive/observable-state.js';
import { classes } from '../_generated/hover-label.module.css.js';
import { HoverLabel } from '../hover-label.js';

// The hover label's placement: in the corner by default (it stays put, nothing follows the pointer), beside the pointer
// when asked.

afterEach(() => {
  document.body.innerHTML = '';
});

function mount(placement?: 'corner' | 'cursor') {
  const area = document.createElement('div');
  document.body.append(area);
  const lines = observableState<readonly string[]>(['PHE 50', '1CBS | A']);
  render(() => <HoverLabel model={{ lines: lines.asReadonly }} area={() => area} placement={placement} />, area);
  return { label: area.querySelector<HTMLElement>('[role=tooltip]')!, area };
}

describe('hover label placement', () => {
  it('is in the corner by default', () => {
    const { label } = mount();
    expect(label.classList.contains(classes.corner)).toBe(true);
  });

  it('follows the pointer only at the cursor', async () => {
    for (const placement of ['corner', 'cursor'] as const) {
      const { label, area } = mount(placement);
      await Promise.resolve();
      area.dispatchEvent(new PointerEvent('pointermove', { clientX: 40, clientY: 30, bubbles: true }));
      await new Promise((r) => setTimeout(r, 40));
      expect(label.classList.contains(classes.corner)).toBe(placement === 'corner');
      // Shown at once (no fade in).
      expect(label.classList.contains(classes.out)).toBe(false);
      if (placement === 'corner') expect(label.style.transform).toBe('');
      document.body.innerHTML = '';
    }
  });
});

describe('hover label content', () => {
  it('renders runs as elements, never as markup', async () => {
    const area = document.createElement('div');
    document.body.append(area);
    const lines = observableState<readonly HoverLine[]>([
      { role: 'title', runs: [{ text: 'PHE', bold: true }, { text: ' <img src=x onerror=alert(1)>' }] },
    ]);
    render(() => <HoverLabel model={{ lines: lines.asReadonly }} area={() => area} />, area);
    const label = area.querySelector<HTMLElement>('[role=tooltip]')!;
    expect(label.querySelector('strong')?.textContent).toBe('PHE');
    expect(label.querySelector('img')).toBeNull();
    expect(label.textContent).toContain('<img src=x onerror=alert(1)>');
  });

  it('keeps the last label a frame when it empties, and swaps in the next without hiding', async () => {
    const area = document.createElement('div');
    document.body.append(area);
    const lines = observableState<readonly HoverLine[]>(['ALA 12']);
    render(() => <HoverLabel model={{ lines: lines.asReadonly }} area={() => area} />, area);
    const label = area.querySelector<HTMLElement>('[role=tooltip]')!;
    lines.set([]);
    expect(label.textContent).toBe('ALA 12');
    // The next label within the same frame: swapped in place.
    lines.set(['GLY 13']);
    expect(label.textContent).toBe('GLY 13');
    lines.set([]);
    await new Promise((r) => setTimeout(r, 100));
    // Hidden (`out`: it fades), its last text kept for the fade.
    expect(label.classList.contains(classes.out)).toBe(true);
    expect(label.textContent).toBe('GLY 13');
  });
});
