/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { render } from '@molstar/prism-ui/dom/render';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { afterEach, describe, expect, it } from 'vitest';
import type { SelectionPart, SelectionSummaryModel } from '../../../model/selection.js';
import type { StyleLayerItem, StyleLayersModel, StyleLayerView } from '../../../model/style-layers.js';
import type { SuperposedEntry } from '../../../model/superposition.js';
import { SelectionSummary } from '../selection-summary.js';
import { StyleLayers } from '../style-layers.js';
import { Superposed } from '../superposed.js';

// The Selection tab's summary and the lists of what its sections made (style layers, superpositions), over mock models.

function mount(view: () => Node) {
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(view, container);
  const buttons = () => [...container.querySelectorAll('button')] as HTMLButtonElement[];
  const button = (text: string) =>
    buttons().find((b) => b.textContent?.trim() === text || b.getAttribute('aria-label') === text);
  return { container, dispose, buttons, button };
}

afterEach(() => {
  document.body.innerHTML = '';
});

describe('SelectionSummary', () => {
  function model() {
    const calls: string[] = [];
    const summary = observableState<string | undefined>(undefined);
    const parts = observableState<readonly SelectionPart[]>([]);
    const selectMode = observableState(false);
    const m: SelectionSummaryModel = {
      summary: summary.asReadonly,
      parts: parts.asReadonly,
      selectMode: selectMode.asReadonly,
      pickLevel: observableState('Residue').asReadonly,
      startSelecting: () => {
        calls.push('select-mode');
        selectMode.set(true);
      },
      focus: () => calls.push('focus'),
      clear: () => calls.push('clear'),
    };
    return { m, calls, summary, parts };
  }

  it('offers Select mode when nothing is selected, then says how to pick', () => {
    const { m, calls } = model();
    const { container, button } = mount(() => <SelectionSummary model={m} />);
    expect(container.textContent).toMatch(/Nothing selected/);
    button('Select on the canvas')!.click();
    expect(calls).toEqual(['select-mode']);
    expect(button('Select on the canvas')).toBeUndefined();
    expect(container.textContent).toMatch(/Click to pick a residue/);
  });

  it('shows the selection, a line per structure (a few, the rest counted), Focus and Clear', () => {
    const { m, calls, summary, parts } = model();
    const { container, button } = mount(() => <SelectionSummary model={m} />);
    summary.set('8 residues');
    parts.set([
      { id: 'a', label: '1CBS', detail: 'chain A · 71 atoms' },
      { id: 'b', label: '1TQN', detail: 'chain A · 30 atoms' },
      { id: 'c', label: '1W0E', detail: 'chain A · 30 atoms' },
      { id: 'd', label: '1W0F', detail: 'chain A · 30 atoms' },
    ]);
    expect(container.textContent).toMatch(/8 residues/);
    expect(container.textContent).toMatch(/1CBSchain A · 71 atoms/);
    expect(container.textContent).not.toMatch(/1W0F/);
    expect(container.textContent).toMatch(/\+1 more structures/);
    button('Focus the selection')!.click();
    button('Clear the selection')!.click();
    expect(calls).toEqual(['focus', 'clear']);
  });
});

describe('StyleLayers', () => {
  function item(key: string, view: StyleLayerView, editable: boolean, log: string[]): StyleLayerItem {
    const form = new ParamsFormModel({ value: PD.Numeric(0.4, { min: 0, max: 1 }) });
    return {
      key,
      view: observableState(view).asReadonly,
      editable,
      editor: () => (editable ? form : undefined),
      remove: () => log.push(`remove ${key}`),
    };
  }

  it('lists the layers with their swatch or kind, edits one inline, removes one, clears all', () => {
    const log: string[] = [];
    const items = observableState<readonly StyleLayerItem[]>([]);
    const m: StyleLayersModel = { items: items.asReadonly, clear: () => log.push('clear') };
    const { container, button, buttons } = mount(() => <StyleLayers model={m} />);
    // Nothing applied: nothing shown.
    expect(container.textContent).toBe('');
    items.set([
      item(
        'c',
        { swatch: '#ff0000', icon: 'paintbrush', where: 'A:10-17', value: '', detail: 'Color #ff0000' },
        true,
        log,
      ),
      item('t', { icon: 'eye-off', where: 'A:200', value: '40% transparent', detail: 'Transparency 0.4' }, true, log),
      item('k', { icon: 'crop', where: 'A:5', value: 'unclipped', detail: 'Excluded' }, false, log),
    ]);
    expect(container.textContent).toMatch(/Applied 3/);
    expect(container.querySelectorAll('[style*="background: #ff0000"]').length).toBe(1);
    // Clipping is not edited in place: two edit buttons for three rows.
    expect(buttons().filter((b) => b.getAttribute('aria-label') === 'Edit')).toHaveLength(2);
    const edit = buttons().filter((b) => b.getAttribute('aria-label') === 'Edit')[1];
    edit.click();
    expect(edit.getAttribute('aria-expanded')).toBe('true');
    expect(container.querySelector('input[type=range]')).not.toBeNull();
    buttons()
      .filter((b) => b.getAttribute('aria-label') === 'Remove')[2]
      .click();
    button('Clear all')!.click();
    expect(log).toEqual(['remove k', 'clear']);
  });
});

describe('Superposed', () => {
  it('lists the superpositions, removes one, and offers Clear all with more than one', () => {
    const log: string[] = [];
    const made = observableState<readonly SuperposedEntry[]>([
      { id: 's1', label: '1W0E', detail: 'onto 1TQN · chains · 0.78 Å', line: '1W0E onto 1TQN by chains, RMSD 0.78 Å' },
    ]);
    const { container, button } = mount(() => (
      <Superposed made={made.asReadonly} remove={(id) => log.push(`remove ${id}`)} clear={() => log.push('clear')} />
    ));
    expect(container.textContent).toMatch(/1W0Eonto 1TQN · chains · 0.78 Å/);
    expect(button('Clear all')).toBeUndefined();
    button('Remove the superposition of 1W0E')!.click();
    made.set([
      ...made.value,
      { id: 's2', label: '1W0F', detail: 'onto 1TQN · chains · 0.78 Å', line: '1W0F onto 1TQN by chains, RMSD 0.78 Å' },
    ]);
    button('Clear all')!.click();
    expect(log).toEqual(['remove s1', 'clear']);
    made.set([]);
    expect(container.textContent).toBe('');
  });
});
