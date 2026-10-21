/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { render } from '@molstar/prism-ui/dom/render';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { SelectionAction, SelectionModel } from '../../../model/selection.js';
import { SelectionBar } from '../selection-bar.js';

// The selection bar: what is selected as the button of a menu of what to do with it, and Clear.

const tick = () => new Promise((r) => setTimeout(r, 0));

function mount() {
  const ran: string[] = [];
  let cleared = 0;
  const summary = observableState<string | undefined>('9 residues');
  const action = (id: string, label: string, hint: string): SelectionAction => ({
    id,
    label,
    hint,
    icon: 'focus',
    run: () => ran.push(id),
  });
  const model: SelectionModel = {
    summary: summary.asReadonly,
    actions: observableState<readonly SelectionAction[]>([
      action('focus', 'Focus', 'Bring the selection into view'),
      action('show-in-tree', 'Show in structure tree', 'Only the rows that hold it'),
      action('more', 'More in the Selection panel…', 'Color, measure, superpose'),
    ]).asReadonly,
    clear: () => cleared++,
  };
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <SelectionBar selection={model} />, container);
  return { container, ran, cleared: () => cleared, summary, dispose };
}

afterEach(() => {
  document.body.innerHTML = '';
});

describe('the selection bar', () => {
  it('opens its menu from the summary: each row a title and a line, a row runs and closes it', async () => {
    const { container, ran, dispose } = mount();
    const button = container.querySelector<HTMLButtonElement>('button[aria-haspopup=menu]')!;
    expect(button.textContent).toContain('9 residues');
    expect(document.querySelector('[role=menu]')).toBeNull();

    button.click();
    await tick();
    expect(button.getAttribute('aria-expanded')).toBe('true');
    const rows = [...document.querySelectorAll<HTMLButtonElement>('[role=menu] [role=menuitem]')];
    expect(rows.map((r) => r.textContent)).toEqual([
      'FocusBring the selection into view',
      'Show in structure treeOnly the rows that hold it',
      'More in the Selection panel…Color, measure, superpose',
    ]);
    rows[1].click();
    await tick();
    expect(ran).toEqual(['show-in-tree']);
    expect(document.querySelector('[role=menu]')).toBeNull();
    dispose();
  });

  it('keeps Clear a button of its own, and is gone with the selection', async () => {
    const { container, cleared, summary, dispose } = mount();
    container.querySelector<HTMLButtonElement>('[aria-label="Clear selection"]')!.click();
    expect(cleared()).toBe(1);
    summary.set(undefined);
    await tick();
    expect(container.querySelector('[role=toolbar]')).toBeNull();
    dispose();
  });

  it('opens on hover and closes a moment after the pointer leaves both the button and the menu', async () => {
    const { container, dispose } = mount();
    const button = container.querySelector<HTMLButtonElement>('button[aria-haspopup=menu]')!;
    const pointer = (type: string, el: Element) =>
      el.dispatchEvent(new PointerEvent(type, { pointerType: 'mouse', bubbles: false }));
    pointer('pointerenter', button);
    await tick();
    expect(document.querySelector('[role=menu]')).not.toBeNull();
    // over to the menu: it stays
    pointer('pointerleave', button);
    const menu = document.querySelector('[role=menu]')!;
    pointer('pointerenter', menu);
    await new Promise((r) => setTimeout(r, 200));
    expect(document.querySelector('[role=menu]')).not.toBeNull();
    // off both: it goes
    pointer('pointerleave', menu);
    await new Promise((r) => setTimeout(r, 200));
    expect(document.querySelector('[role=menu]')).toBeNull();
    dispose();
  });

  it('runs a row pressed inside a shadow root (a press outside the bar closes it, one inside does not)', async () => {
    const host = document.createElement('div');
    document.body.append(host);
    const shadow = host.attachShadow({ mode: 'open' });
    const container = document.createElement('div');
    shadow.append(container);
    const ran: string[] = [];
    const model: SelectionModel = {
      summary: observableState<string | undefined>('REA 200').asReadonly,
      actions: observableState<readonly SelectionAction[]>([
        { id: 'focus', label: 'Focus', icon: 'focus', run: () => ran.push('focus') },
      ]).asReadonly,
      clear: () => {},
    };
    const dispose = render(() => <SelectionBar selection={model} />, container);
    shadow.querySelector<HTMLButtonElement>('button[aria-haspopup=menu]')!.click();
    await tick();
    const row = shadow.querySelector<HTMLButtonElement>('[role=menuitem]')!;
    row.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, composed: true }));
    await tick();
    // still open: the press was inside the bar's menu, though the document sees the shadow host as its target
    expect(shadow.querySelector('[role=menu]')).not.toBeNull();
    row.click();
    expect(ran).toEqual(['focus']);
    // a press outside closes it
    shadow.querySelector<HTMLButtonElement>('button[aria-haspopup=menu]')!.click();
    await tick();
    document.body.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, composed: true }));
    await tick();
    expect(shadow.querySelector('[role=menu]')).toBeNull();
    dispose();
  });
});
