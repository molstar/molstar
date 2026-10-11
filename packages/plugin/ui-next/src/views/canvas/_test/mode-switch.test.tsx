/** @jest-environment jsdom */
/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { type CanvasMode, ModeSwitch } from '../mode-switch.js';

describe('ModeSwitch', () => {
  it('shows the mode, switches it, and offers granularity only while selecting', () => {
    const mode = observableState<CanvasMode>('interact');
    const granularity = observableState('residue');
    const container = document.createElement('div');
    render(
      () => (
        <ModeSwitch
          mode={mode.asReadonly}
          setMode={(m) => mode.set(m)}
          granularity={granularity.asReadonly}
          granularities={[
            ['residue', 'Residue'],
            ['chain', 'Chain'],
          ]}
          setGranularity={(g) => granularity.set(g)}
        />
      ),
      container,
    );
    const [interact, select] = container.querySelectorAll('button');
    expect(interact.getAttribute('aria-pressed')).toBe('true');
    expect(container.querySelector('select')).toBeNull();

    select.click();
    expect(mode.value).toBe('select');
    expect(select.getAttribute('aria-pressed')).toBe('true');
    const by = container.querySelector('select') as HTMLSelectElement;
    expect(by.options[by.selectedIndex].text).toBe('Residue');
    by.value = '1';
    by.dispatchEvent(new Event('change'));
    expect(granularity.value).toBe('chain');

    interact.click();
    expect(container.querySelector('select')).toBeNull();
  });
});
