// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { createRoot } from '../../../dom/scope.js';
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
    // In a scope, as the shell installs them (its listeners are scoped).
    uninstall = createRoot((dispose) => {
      const off = installTooltips(root);
      return () => {
        off();
        dispose();
      };
    });
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

  it("draws a cut-off tree row's card at the row's height, each part in its own size and weight", () => {
    root.insertAdjacentHTML(
      'beforeend',
      '<div id="row" data-tip-row=""><span data-tip-part="" style="font-size: 12px; font-weight: 500">PAE of a model</span>' +
        '<span data-tip-part="" style="font-size: 11px">1098 residues</span></div>',
    );
    const row = root.querySelector('#row') as HTMLElement;
    const label = row.firstElementChild as HTMLElement;
    // jsdom has no layout: the label is cut off, the row 26 px high, the card's line 16.8 px.
    Object.defineProperty(label, 'scrollWidth', { value: 120 });
    Object.defineProperty(label, 'clientWidth', { value: 60 });
    row.getBoundingClientRect = () =>
      ({ top: 100, height: 26, left: 0, right: 200, bottom: 126, width: 200 }) as DOMRect;
    tip().style.lineHeight = '16.8px';
    over(label);
    vi.advanceTimersByTime(400);
    expect(tip().hidden).toBe(false);
    const [a, b] = [...tip().querySelectorAll('span')];
    expect([a.textContent, a.style.fontSize, a.style.fontWeight]).toEqual(['PAE of a model', '12px', '500']);
    expect([b.textContent, b.style.fontSize]).toEqual(['1098 residues', '11px']);
    expect(Number.parseFloat(tip().style.paddingBlock)).toBeCloseTo(4.6);
    // Wrapped, taller than the row: the card keeps its own padding.
    label.dispatchEvent(new PointerEvent('pointerout', { bubbles: true }));
    tip().getBoundingClientRect = () =>
      ({ top: 0, height: 60, left: 0, right: 300, bottom: 60, width: 300 }) as DOMRect;
    over(label);
    vi.advanceTimersByTime(400);
    expect(tip().hidden).toBe(false);
    expect(tip().style.paddingBlock).toBe('');
  });
});
