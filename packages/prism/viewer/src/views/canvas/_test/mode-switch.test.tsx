/** @jest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { type CanvasMode, ModeSwitch, type ModeSwitchProps, PickLevel } from '../mode-switch.js';

describe('ModeSwitch', () => {
  it('shows the mode, switches it, and offers granularity only while selecting', () => {
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
    const [interact, select] = container.querySelectorAll('button');
    const level = () => container.querySelector('[aria-haspopup]') as HTMLButtonElement | null;
    expect(interact.getAttribute('aria-pressed')).toBe('true');
    expect(level()).toBeNull();

    // Selecting shows the level below the switch; the switch's own buttons stay as they were.
    select.click();
    expect(mode.value).toBe('select');
    expect(select.getAttribute('aria-pressed')).toBe('true');
    expect(level()!.textContent).toContain('Residue');
    level()!.click();
    const rows = [...container.querySelectorAll('[role=menuitemradio]')] as HTMLButtonElement[];
    expect(rows.map((r) => r.getAttribute('aria-checked'))).toEqual(['true', 'false']);
    rows[1].click();
    expect(granularity.value).toBe('chain');
    expect(container.querySelector('[role=menu]')).toBeNull();
    expect(level()!.textContent).toContain('Chain');

    // Escape closes the list.
    level()!.click();
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
    expect(container.querySelector('[role=menu]')).toBeNull();

    interact.click();
    expect(level()).toBeNull();
    container.remove();
  });
});
