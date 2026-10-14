// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { installTooltips } from '../tooltip.js';

// A tip whose anchor goes (a panel closed under a resting pointer) must not stay on screen.

describe('tooltips', () => {
  let root: HTMLElement;
  let uninstall: () => void;
  beforeEach(() => {
    vi.useFakeTimers();
    root = document.createElement('div');
    document.body.append(root);
    root.innerHTML = '<div id="panel"><button id="a" data-tip="Reset">x</button></div>';
    uninstall = installTooltips(root);
  });
  afterEach(() => {
    uninstall();
    root.remove();
    vi.useRealTimers();
  });
  const tip = () => root.querySelector<HTMLElement>('[role="tooltip"]') as HTMLElement;
  const over = (el: Element) => el.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
  const a = () => root.querySelector('#a') as HTMLElement;

  it('shows after the delay and hides on pointerout', () => {
    over(a());
    expect(tip().hidden).toBe(true);
    vi.advanceTimersByTime(400);
    expect(tip().hidden).toBe(false);
    expect(tip().textContent).toBe('Reset');
    a().dispatchEvent(new PointerEvent('pointerout', { bubbles: true }));
    expect(tip().hidden).toBe(true);
  });

  it('hides when its anchor is removed, at the next pointer event in the root', () => {
    over(a());
    vi.advanceTimersByTime(400);
    expect(tip().hidden).toBe(false);
    a().remove();
    root.dispatchEvent(new PointerEvent('pointermove', { bubbles: true }));
    expect(tip().hidden).toBe(true);
  });

  it('is not shown for an anchor removed while the pointer rested', () => {
    over(a());
    a().remove();
    vi.advanceTimersByTime(400);
    expect(tip().hidden).toBe(true);
  });

  it('is never shown for an anchor in an inert part (a hidden panel)', () => {
    document.getElementById('panel')!.setAttribute('inert', '');
    over(a());
    vi.advanceTimersByTime(400);
    expect(tip().hidden).toBe(true);
  });

  it('ends when the pointer moves over something that is not the anchor', () => {
    over(a());
    vi.advanceTimersByTime(400);
    expect(tip().hidden).toBe(false);
    root.dispatchEvent(new PointerEvent('pointermove', { bubbles: true }));
    expect(tip().hidden).toBe(true);
  });
});
