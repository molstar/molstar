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
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { CommandManager } from '../../../commands/manager.js';
import { render } from '../../../dom/render.js';
import { flush, RuntimeProbe } from '../../../dom/scope.js';
import { observableState } from '../../../reactive/observable-state.js';
import { CommandMenu } from '../command-menu.js';

// The form mode as a view of the invocation's state (rendering.md §2.6): the input, the run button, the field tokens
// and the options are each written by a `watch`, which listens to what decides its output.

let rendered = 0;
beforeEach(() => {
  rendered = 0;
  RuntimeProbe.on = (kind) => {
    if (kind === 'render') rendered++;
  };
});
afterEach(() => {
  RuntimeProbe.on = undefined;
  document.body.innerHTML = '';
});

function open() {
  const manager = new CommandManager<object>({
    context: observableState<object>({}).asReadonly,
    formats: { get: () => undefined },
  });
  manager.register([
    {
      id: 'look.show',
      title: 'Show as',
      category: 'Look',
      params: () => ({
        type: PD.Select('cartoon', PD.arrayToOptions(['cartoon', 'surface', 'spacefill'])),
        size: PD.Numeric(1, { min: 0, max: 5, step: 1 }),
      }),
      run: () => ({}),
    },
  ]);
  manager.open('look.show');
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <CommandMenu manager={manager} onClose={vi.fn()} />, container);
  const input = () => container.querySelector<HTMLInputElement>('input[role=combobox]')!;
  const options = () => [...container.querySelectorAll<HTMLElement>('[role=option]')];
  const key = (k: string) => {
    input().dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true }));
    flush();
  };
  return { manager, container, input, options, key };
}

describe('the command menu’s form mode', () => {
  it('the options and the input follow the pick, the ids stay', () => {
    const t = open();
    flush();
    const ids = t.options().map((o) => o.id);
    expect(ids.length).toBe(3);
    expect(t.input().getAttribute('aria-expanded')).toBe('true');
    const picked = () => t.options().map((o) => o.getAttribute('aria-selected'));
    // Nothing picked yet: no active descendant.
    expect(picked()).toEqual(['false', 'false', 'false']);
    expect(t.input().hasAttribute('aria-activedescendant')).toBe(false);

    // A move of the pick: the options' marks and the input's `aria-activedescendant` follow, ids stay.
    t.key('ArrowDown');
    const now = picked().indexOf('true');
    expect(now).toBe(0);
    expect(t.input().getAttribute('aria-activedescendant')).toBe(ids[now]);
    expect(t.options().map((o) => o.id)).toEqual(ids);
    expect(picked().filter((x) => x === 'true')).toHaveLength(1);
  });

  it('the sentence’s tokens name the active field and the values; the input names the field', () => {
    const t = open();
    flush();
    const tokens = () => [...t.container.querySelectorAll<HTMLElement>('button[data-empty], button[class*=value]')];
    const active = () => t.container.querySelector('[data-active]')?.textContent?.trim();
    expect(t.input().getAttribute('aria-label')).toBeTruthy();
    const field = t.input().getAttribute('aria-label');
    expect(active()).toBe(field);
    expect(tokens().length).toBeGreaterThan(0);

    // Typing into the field writes the input (not the other way): a value that normalizes comes back as the state's.
    t.input().value = 'surf';
    t.input().dispatchEvent(new Event('input', { bubbles: true }));
    flush();
    expect(t.input().value).toBe('surf');
    expect(t.options().some((o) => o.textContent?.toLowerCase().includes('surface'))).toBe(true);
  });

  it('a keystroke renders the input and the options’ watches that concern it, no more', () => {
    const t = open();
    flush();
    const before = rendered;
    t.key('ArrowDown');
    // The input; the two options whose mark changed (the pick is their source, not the whole state).
    expect(rendered - before).toBeLessThanOrEqual(1 + t.options().length);
  });
});
