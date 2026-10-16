/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { render } from '@molstar/prism-ui/dom/render';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { type CanvasMode, ModeSwitch, type ModeSwitchProps, PickLevel } from '../mode-switch.js';

describe('ModeSwitch', () => {
  it('shows the mode and switches it; the picking levels are listed directly, one click from Select', () => {
    const mode = observableState<CanvasMode>('interact');
    const granularity = observableState('residue');
    const container = document.createElement('div');
    const props: ModeSwitchProps = {
      mode: mode.asReadonly,
      setMode: (m) => mode.set(m),
      granularity: granularity.asReadonly,
      granularities: [
        { value: 'residue', label: 'Residue', icon: 'pick-residue' },
        { value: 'chain', label: 'Chain', icon: 'pick-chain' },
        { value: 'chainInstances', label: 'Chain', icon: 'pick-chain', group: 'All instances' },
      ],
      setGranularity: (g) => granularity.set(g),
    };
    render(
      () => (
        <>
          <ModeSwitch {...props} />
          <PickLevel {...props} />
        </>
      ),
      container,
    );
    document.body.append(container);
    const [interact, select] = container.querySelectorAll<HTMLButtonElement>('.msp-viewer-mode-switch-segment');
    const rows = () => [...container.querySelectorAll('[role=menuitemradio]')] as HTMLButtonElement[];
    expect(interact.getAttribute('aria-pressed')).toBe('true');

    // The list is there (shown by CSS on hover or focus of Select): no button to open first, the current one marked;
    // Select controls it (the keyboard's way in).
    expect(select.getAttribute('aria-controls')).toBe(container.querySelector('[role=menu]')?.id);
    expect(rows().map((r) => r.getAttribute('aria-checked'))).toEqual(['true', 'false', 'false']);
    expect(container.querySelector('[role=menu]')?.textContent).toContain('All instances');

    // Picking a level sets it and turns select mode on (it was off).
    rows()[1].click();
    expect(granularity.value).toBe('chain');
    expect(mode.value).toBe('select');
    expect(select.getAttribute('aria-pressed')).toBe('true');
    expect(rows().map((r) => r.getAttribute('aria-checked'))).toEqual(['false', 'true', 'false']);

    // Another pick keeps select mode on; Interact turns it off and keeps the level.
    rows()[2].click();
    expect(granularity.value).toBe('chainInstances');
    expect(mode.value).toBe('select');
    interact.click();
    expect(mode.value).toBe('interact');
    expect(granularity.value).toBe('chainInstances');
    container.remove();
  });

  it('is reached by the keyboard: ↑ on Select goes to the current level, the arrows move, Escape goes back', () => {
    const mode = observableState<CanvasMode>('interact');
    const granularity = observableState('chain');
    const container = document.createElement('div');
    const props: ModeSwitchProps = {
      mode: mode.asReadonly,
      setMode: (m) => mode.set(m),
      granularity: granularity.asReadonly,
      granularities: [
        { value: 'residue', label: 'Residue', icon: 'pick-residue' },
        { value: 'chain', label: 'Chain', icon: 'pick-chain' },
        { value: 'element', label: 'Atom', icon: 'pick-residue' },
      ],
      setGranularity: (g) => granularity.set(g),
    };
    render(
      () => (
        <>
          <ModeSwitch {...props} />
          <PickLevel {...props} />
        </>
      ),
      container,
    );
    document.body.append(container);
    const select = container.querySelector<HTMLButtonElement>('[aria-label="Select"]')!;
    const rows = [...container.querySelectorAll<HTMLButtonElement>('[role=menuitemradio]')];
    const key = (el: Element, k: string) =>
      el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true }));
    select.focus();
    key(select, 'ArrowUp');
    expect(document.activeElement).toBe(rows[1]);
    expect(select.getAttribute('aria-expanded')).toBe('true');
    key(rows[1], 'ArrowDown');
    expect(document.activeElement).toBe(rows[2]);
    key(rows[2], 'ArrowDown');
    expect(document.activeElement).toBe(rows[0]);
    key(rows[0], 'End');
    expect(document.activeElement).toBe(rows[2]);
    key(rows[2], 'Escape');
    expect(document.activeElement).toBe(select);
    expect(select.getAttribute('aria-expanded')).toBe('false');
    container.remove();
  });

  it('lets the focus go on a tap outside, so the list closes on a touch screen (where the focus shows it)', () => {
    const mode = observableState<CanvasMode>('interact');
    const granularity = observableState('residue');
    const container = document.createElement('div');
    const props: ModeSwitchProps = {
      mode: mode.asReadonly,
      setMode: (m) => mode.set(m),
      granularity: granularity.asReadonly,
      granularities: [{ value: 'residue', label: 'Residue', icon: 'pick-residue' }],
      setGranularity: (g) => granularity.set(g),
    };
    render(
      () => (
        <>
          <ModeSwitch {...props} />
          <PickLevel {...props} />
        </>
      ),
      container,
    );
    const canvas = document.createElement('canvas');
    document.body.append(container, canvas);
    const select = container.querySelector<HTMLButtonElement>('[aria-controls]')!;
    select.focus();
    expect(document.activeElement).toBe(select);
    // A tap on Select itself keeps it.
    select.dispatchEvent(new Event('pointerdown', { bubbles: true, composed: true }));
    expect(document.activeElement).toBe(select);
    // A tap on the canvas (not focusable) lets it go.
    canvas.dispatchEvent(new Event('pointerdown', { bubbles: true, composed: true }));
    expect(document.activeElement).not.toBe(select);
    container.remove();
    canvas.remove();
  });
});
