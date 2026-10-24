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

// The console's view over its model (_spec §9): Enter runs the input, ⇧Enter does not, ↑ recalls, Tab completes. The
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
    expect(t.input.value).toBe('show sticks');
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
    // After a space, the next word's completions take the line (one row: the composer keeps its height).
    t.type('show ');
    await tick();
    const chip = [...t.container.querySelectorAll('button')].find((b) => b.textContent === 'sticks')!;
    chip.click();
    expect(t.input.value).toBe('show sticks');
  });
});
