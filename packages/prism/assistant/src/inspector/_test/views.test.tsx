/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { render } from '@molstar/prism-ui/dom/render';
import { flush } from '@molstar/prism-ui/dom/scope';
import { mockInspector } from '../../testing/mock-inspector.js';
import { findLimit, type InspectorModel } from '../model.js';
import { InspectorBody, InspectorHeader } from '../page.js';

function mount(model: InspectorModel) {
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(
    () => (
      <>
        <header>
          <InspectorHeader model={model} />
        </header>
        <InspectorBody model={model} />
      </>
    ),
    container,
  );
  const selected = () => container.querySelector<HTMLElement>('[role=listbox] [aria-selected=true]')?.dataset.row;
  // A key moves the selection; the rows show it after the render the change queued (rendering.md §2.3).
  const key = (key: string) => {
    document.body.dispatchEvent(new KeyboardEvent('keydown', { key, bubbles: true }));
    flush();
  };
  return { container, dispose, selected, key };
}

afterEach(() => {
  document.body.innerHTML = '';
});

describe('Inspector views', () => {
  it('shows the trace: turns, requests and steps as rows; j/k and the arrows move the selection', () => {
    const model = mockInspector();
    const { container, selected, key } = mount(model);
    const rows = [...container.querySelectorAll<HTMLElement>('[data-row]')].map((r) => r.dataset.row);
    expect(rows).toContain('turn:turn-1');
    expect(rows).toContain('req:turn-2');
    expect(rows.filter((r) => r?.startsWith('step:')).length).toBeGreaterThan(20);
    // The failed tool call and the nested helper are marked.
    expect(container.querySelector('[role=listbox] [data-failed]')).not.toBeNull();
    expect(container.querySelector('[role=listbox] [data-nested]')).not.toBeNull();
    // The overview strip: a block per step, the selected one pressed.
    expect(container.querySelectorAll('[aria-pressed=true]')).toHaveLength(1);

    const last = rows.at(-1);
    expect(selected()).toBe(last);
    key('k');
    expect(selected()).toBe(rows.at(-2));
    key('ArrowUp');
    expect(selected()).toBe(rows.at(-3));
    key('j');
    expect(selected()).toBe(rows.at(-2));
    // The tab counts: 3 turns, 14 tools.
    const tabs = [...container.querySelectorAll('header [role=tab]')].map((t) => t.textContent);
    expect(tabs).toEqual(['Trace3', 'Tools14', 'Limits']);
  });

  it('shows a model call’s window and switches detail tabs; a tool call’s definition leads to Tools', () => {
    const model = mockInspector();
    const { container } = mount(model);
    model.focus('step:s11');
    expect(container.textContent).toContain('It does not extend the previous request');
    const tab = (label: string) =>
      [...container.querySelectorAll<HTMLElement>('[role=tab]')].find((t) => t.textContent === label)!;
    tab('Timing').click();
    expect(container.textContent).toContain('until the first output');
    model.focus('step:s14');
    expect(container.textContent).toContain('Could not read the expression');
    tab('Definition').click();
    [...container.querySelectorAll<HTMLElement>('button')].find((b) => b.textContent?.includes('in Tools'))!.click();
    expect(model.view.value).toBe('tools');
    expect(container.querySelector('[role=listbox] [aria-selected=true]')?.getAttribute('data-row')).toBe(
      'tool:selection.set',
    );
  });

  it('lists the tools by group and moves through them; a tool’s input as fields', () => {
    const model = mockInspector();
    model.focus('tools');
    const { container, selected, key } = mount(model);
    expect(selected()).toBe('group:core');
    expect(container.querySelector('[data-row="group:code"]')?.textContent).toContain('unavailable here');
    key('j');
    expect(selected()).toBe('tool:scene.get');
    model.openTool('structure.load');
    expect(container.textContent).toContain('asks the person first');
    expect(container.textContent).toContain('One of: pdb, afdb, url');
  });

  it('commits a limit on Enter, puts it back on Escape, and resets it to the default', () => {
    const model = mockInspector();
    model.focus('limits');
    const { container } = mount(model);
    const input = container.querySelector<HTMLInputElement>('input[aria-label="Model calls"]')!;
    expect(input.value).toBe('24');
    const press = (key: string) => input.dispatchEvent(new KeyboardEvent('keydown', { key, bubbles: true }));

    input.value = '40';
    press('Enter');
    expect(findLimit(model.limits.value, 'modelCalls')?.value).toBe(40);
    expect(input.value).toBe('40');

    input.value = '33';
    press('Escape');
    expect(input.value).toBe('40');
    expect(findLimit(model.limits.value, 'modelCalls')?.value).toBe(40);

    // Over the ceiling: refused, the field shows the value as it is.
    input.value = '5k';
    press('Enter');
    expect(input.value).toBe('40');

    // The way back is there while the value differs from the default.
    const shown = () =>
      [...container.querySelectorAll('button')].some((b) => b.textContent === 'Reset' && b.style.display !== 'none');
    flush();
    expect(shown()).toBe(true);
    const reset = [...container.querySelectorAll<HTMLButtonElement>('button')].find((b) => b.textContent === 'Reset')!;
    reset.click();
    expect(findLimit(model.limits.value, 'modelCalls')?.value).toBe(24);
    flush();
    expect(input.value).toBe('24');
    expect(shown()).toBe(false);
  });

  it('says so before the agent has run', () => {
    const { container } = mount(mockInspector({ empty: true }));
    expect(container.textContent).toContain('The agent has not run yet');
  });
});
