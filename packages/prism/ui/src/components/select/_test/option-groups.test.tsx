/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { OptionList, type SelectOption } from '../select.js';

// Options with groups (Mol*'s categories): the ungrouped first, then each group a fold, the current option's open.

afterEach(() => {
  document.body.innerHTML = '';
});

const options: SelectOption[] = [
  ['uniform', 'Uniform'],
  ['element-symbol', 'Element Symbol', 'Atom Property'],
  ['chain-id', 'Chain Id', 'Chain Property'],
  ['entity-id', 'Entity Id', 'Chain Property'],
];

describe('an option list with groups', () => {
  it('folds the groups, the current option’s open; a head opens its group', async () => {
    const container = document.createElement('div');
    document.body.append(container);
    const pick = vi.fn();
    render(
      () => (
        <OptionList options={options} value={observableState('chain-id').asReadonly} onPick={pick} onClose={() => {}} />
      ),
      container,
    );
    const heads = [...container.querySelectorAll<HTMLButtonElement>('[data-group-head]')];
    expect(heads.map((h) => h.textContent)).toEqual(['Atom Property', 'Chain Property']);
    expect(heads.map((h) => h.getAttribute('aria-expanded'))).toEqual(['false', 'true']);
    const shown = () =>
      [...container.querySelectorAll('[role="option"]')]
        .filter((o) => !o.closest('[hidden]'))
        .map((o) => o.textContent);
    expect(shown()).toEqual(['Uniform', 'Chain Id', 'Entity Id']);
    heads[0].click();
    expect(shown()).toEqual(['Uniform', 'Element Symbol', 'Chain Id', 'Entity Id']);
    container.querySelectorAll<HTMLButtonElement>('[role="option"]')[1].click();
    expect(pick).toHaveBeenCalledWith('element-symbol');
  });
});
