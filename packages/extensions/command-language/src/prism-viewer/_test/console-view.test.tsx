/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui';
import type { ViewerApi } from '@molstar/prism-viewer/api/area';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { EmptyEnv } from '../../language/env.js';
import { ConsoleView } from '../console.js';
import { ConsoleModel } from '../console-model.js';

// The console's view over its model (_spec §9): Enter runs the input, ⇧Enter does not, ↑ recalls, Tab completes, the
// completion list takes ↑ ↓ Enter Esc while open. The
// API is a stand-in: what the view does with a run's result, not the run.

afterEach(() => {
  document.body.innerHTML = '';
});

const tick = () => new Promise((r) => setTimeout(r, 0));

function open() {
  const run = vi.fn(async () => ({
    ok: true,
    did: 'Did it.',
    changed: {},
    value: { outputs: [{ kind: 'value', text: '42' }], ran: 1 },
  }));
  const api = {
    script: {
      env: () => EmptyEnv,
      preview: () => ({
        dialect: 'pymol',
        commands: [{ kind: 'read', read: { kind: 'names' } }],
        diagnostics: [],
        session: {},
      }),
      run,
    },
    history: {},
  } as unknown as ViewerApi;
  const model = new ConsoleModel(api);
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <ConsoleView model={model} />, container);
  const input = container.querySelector('textarea')!;
  const key = (k: string, shift = false) =>
    input.dispatchEvent(new KeyboardEvent('keydown', { key: k, shiftKey: shift, bubbles: true, cancelable: true }));
  const type = (text: string) => {
    input.value = text;
    input.dispatchEvent(new Event('input', { bubbles: true }));
  };
  return { container, input, key, type, run, model };
}

describe('the console view', () => {
  it('runs on Enter, shows what it printed, and recalls with ↑', async () => {
    const t = open();
    t.type('count_atoms');
    t.key('Enter', true);
    expect(t.run).not.toHaveBeenCalled();
    t.key('Enter');
    await tick();
    await tick();
    expect(t.run).toHaveBeenCalledWith({ dialect: 'pymol', text: 'count_atoms', files: [] });
    expect(t.container.textContent).toContain('42');
    expect(t.input.value).toBe('');
    t.key('ArrowUp');
    expect(t.input.value).toBe('count_atoms');
  });

  it('completes the input with Tab', () => {
    const t = open();
    t.type('show stic');
    t.key('Tab');
    // A PyMOL argument with more after it comes with its comma.
    expect(t.input.value).toBe('show sticks, ');
  });

  it('lists the completions: ↓ picks, Enter takes the one picked, Esc closes, Enter with none picked runs', async () => {
    const t = open();
    t.type('show ');
    const list = () => t.container.querySelector('[role=listbox]');
    const options = () => [...t.container.querySelectorAll('[role=option]')];
    expect(list()).not.toBeNull();
    expect(t.input.getAttribute('aria-controls')).toBe(list()!.id);
    // The signature says which argument is being typed.
    expect(t.container.querySelector('strong')?.textContent).toBe('[representation]');
    t.key('ArrowDown');
    t.key('ArrowDown');
    expect(options()[1].getAttribute('aria-selected')).toBe('true');
    expect(t.input.getAttribute('aria-activedescendant')).toBe(options()[1].id);
    const second = options()[1].querySelector('[class*=label]')!.textContent!;
    t.key('Enter');
    expect(t.run).not.toHaveBeenCalled();
    expect(t.input.value).toBe(`show ${second}, `);
    t.key('Escape');
    expect(list()).toBeNull();
    expect(t.input.hasAttribute('aria-controls')).toBe(false);
    t.type('count_atoms');
    t.key('Enter');
    await tick();
    expect(t.run).toHaveBeenCalledWith({ dialect: 'pymol', text: 'count_atoms', files: [] });
  });

  it('lists every command on Ctrl+I in the empty input, until the next keystroke; Tab there moves the focus on', () => {
    const t = open();
    const options = () => t.container.querySelectorAll('[role=option]');
    const listKey = () =>
      t.input.dispatchEvent(
        new KeyboardEvent('keydown', { key: 'i', code: 'KeyI', ctrlKey: true, bubbles: true, cancelable: true }),
      );
    // Tab is the page's: not taken.
    expect(t.input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Tab', bubbles: true, cancelable: true }))).toBe(
      true,
    );
    expect(options()).toHaveLength(0);
    listKey();
    expect(options().length).toBeGreaterThan(12);
    t.key('Escape');
    expect(options()).toHaveLength(0);
    listKey();
    t.key('Tab');
    expect(t.input.value).toMatch(/^\w+ ?$/);
    t.type('');
    expect(options()).toHaveLength(0);
  });

  it('picks none when the rows change under a pick (the dialect switched): Enter runs, nothing throws', async () => {
    const t = open();
    t.type('color ');
    t.key('ArrowUp');
    expect(t.input.hasAttribute('aria-activedescendant')).toBe(true);
    t.model.setDialect('vmd');
    expect(t.input.hasAttribute('aria-activedescendant')).toBe(false);
    t.key('Enter');
    await tick();
    expect(t.run).toHaveBeenCalled();
  });

  it('goes on through the history when a recalled line has completions', async () => {
    const t = open();
    for (const line of ['show sticks', 'color ']) {
      t.type(line);
      t.key('Escape');
      t.key('Enter');
      await tick();
      await tick();
    }
    t.key('ArrowUp');
    expect(t.input.value).toBe('color ');
    expect(t.container.querySelector('[role=listbox]')).toBeNull();
    t.key('ArrowUp');
    expect(t.input.value).toBe('show sticks');
  });

  it("shows the command's help over the composer while it is typed, without a list where nothing completes", () => {
    const t = open();
    const card = () => t.container.querySelector('[class*=card]');
    t.type('fetch 1abc, ');
    expect(card()?.textContent).toContain('Load an entry');
    expect(t.container.querySelector('[role=listbox]')).toBeNull();
    expect(card()?.querySelector('strong')?.textContent).toBe('[name]');
    t.key('Escape');
    expect(card()).toBeNull();
    t.type('');
    expect(card()).toBeNull();
  });

  it('marks the part typed and shows a color as its swatch', () => {
    const t = open();
    t.type('color sal');
    const option = t.container.querySelector('[role=option]')!;
    expect(option.querySelector('mark')?.textContent).toBe('sal');
    expect(option.querySelector<HTMLElement>('[class*=swatch]')?.getAttribute('style')).toMatch(/background: #/);
  });

  it('says what the line will do, offers its completions, and fills an example', async () => {
    const t = open();
    expect(t.container.textContent).toContain('PyMOL and VMD commands');
    const example = [...t.container.querySelectorAll('button')].find((b) => b.textContent === 'show sticks, resn REA')!;
    example.click();
    expect(t.input.value).toBe('show sticks, resn REA');
    t.type('get_names');
    await tick();
    expect(t.container.textContent).toContain('PyMOL · List the names');
    // After a space, the next word's completions list over the composer; a click takes one, the input keeps focus.
    t.type('show ');
    await tick();
    const option = [...t.container.querySelectorAll('[role=option]')].find((o) => o.textContent?.includes('sticks'))!;
    option.dispatchEvent(new MouseEvent('mousedown', { bubbles: true, cancelable: true }));
    expect(t.input.value).toBe('show sticks, ');
  });
});
