// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { createRoot, runInScope, ViewScope } from '../../../dom/index.js';
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

  it("puts a row's card outside its edge: on the side the edge asks for, else the roomier side", () => {
    const rect = (left: number, right: number, top = 100, height = 24) =>
      ({ left, right, top, bottom: top + height, width: right - left, height, x: left, y: top }) as DOMRect;
    root.insertAdjacentHTML(
      'beforeend',
      '<div id="edge" data-tip-edge=""><div id="r" data-tip-row=""><span data-tip-part="">A long name</span></div></div>',
    );
    const edge = root.querySelector('#edge') as HTMLElement;
    const row = root.querySelector('#r') as HTMLElement;
    const part = row.firstElementChild as HTMLElement;
    Object.defineProperty(part, 'scrollWidth', { value: 120 });
    Object.defineProperty(part, 'clientWidth', { value: 60 });
    vi.stubGlobal('innerWidth', 1000);
    row.getBoundingClientRect = () => rect(700, 940);
    tip().getBoundingClientRect = () => rect(0, 200, 0, 24);
    const x = () => Number(/translate\((-?[\d.]+)px/.exec(tip().style.transform)?.[1]);
    const showAt = (left: number, right: number, ask = '') => {
      part.dispatchEvent(new PointerEvent('pointerout', { bubbles: true }));
      edge.dataset.tipEdge = ask;
      edge.getBoundingClientRect = () => rect(left, right, 0, 600);
      over(part);
      vi.advanceTimersByTime(400);
    };
    // A panel at the right: the room is on its left (toward the main area).
    showAt(680, 960);
    expect(tip().hidden).toBe(false);
    expect(x()).toBe(680 - 8 - 200);
    // A panel at the left: its right.
    showAt(0, 300);
    expect(x()).toBe(308);
    // Asked for the right, with room there: the right, though the left has more.
    showAt(500, 700, 'right');
    expect(x()).toBe(708);
  });

  it("shows a row's card when the keyboard focuses the row (a tree's arrows), not a button in it", () => {
    root.insertAdjacentHTML(
      'beforeend',
      '<div id="k" data-tip-row="" tabindex="0"><span data-tip-part="">A long name</span><button id="kb">x</button></div>',
    );
    const row = root.querySelector('#k') as HTMLElement;
    const part = row.firstElementChild as HTMLElement;
    Object.defineProperty(part, 'scrollWidth', { value: 120 });
    Object.defineProperty(part, 'clientWidth', { value: 60 });
    // jsdom has no :focus-visible; the row matches it here.
    const matches = row.matches.bind(row);
    row.matches = ((s: string) => (s === ':focus-visible' ? true : matches(s))) as typeof row.matches;
    row.focus();
    expect(tip().hidden).toBe(false);
    (root.querySelector('#kb') as HTMLElement).focus();
    expect(tip().hidden).toBe(true);
  });

  it('waits for its delay on the scope: none shows while the scope is paused, the wait resumes with it', () => {
    const scoped = document.createElement('div');
    document.body.append(scoped);
    scoped.innerHTML = '<button id="b" data-tip="Center">x</button>';
    const scope = new ViewScope();
    const off = runInScope(scope, () => installTooltips(scoped));
    const b = scoped.querySelector('#b')!;
    const tipB = () => scoped.querySelector<HTMLElement>('[role="tooltip"]')!;
    over(b);
    scope.setActive(false);
    vi.advanceTimersByTime(400);
    expect(tipB().hidden).toBe(true);
    scope.setActive(true);
    vi.advanceTimersByTime(400);
    expect(tipB().hidden).toBe(false);
    off();
    scope.dispose();
    scoped.remove();
  });
});
