/** @jest-environment jsdom */
/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mockCommands } from '../../../_stories/mocks/commands.js';
import { render } from '../../../dom/render.js';
import { ShellModel } from '../../../shell/shell.js';
import { CommandMenu } from '../command-menu.js';

const flush = () => new Promise((r) => setTimeout(r, 0));

function mount() {
  const shell = new ShellModel();
  const manager = mockCommands(shell);
  const closed = jest.fn();
  const results: string[] = [];
  manager.outcomes.subscribe((o) => {
    if (o.kind === 'done') results.push(o.result.message ?? '');
  });
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <CommandMenu manager={manager} onClose={closed} />, container);
  const input = () => container.querySelector('input') as HTMLInputElement;
  const type = (text: string) => {
    input().value = text;
    input().dispatchEvent(new Event('input'));
  };
  const key = (k: string) => input().dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  const rows = () => [...container.querySelectorAll('[role=option]')].map((r) => r.textContent ?? '');
  return { manager, closed, results, container, input, type, key, rows };
}

afterEach(() => {
  document.body.innerHTML = '';
});

describe('CommandMenu', () => {
  it('lists commands, searches, and runs one without arguments', async () => {
    const p = mount();
    expect(p.rows().length).toBeGreaterThan(10);
    p.type('reset zo');
    expect(p.rows()[0]).toContain('Reset zoom');
    p.key('Enter');
    await flush();
    expect(p.results).toEqual(['Zoom reset']);
    expect(p.closed).toHaveBeenCalled();
  });

  it('reads an alias with arguments into a sentence, completes, and runs', async () => {
    const p = mount();
    p.type('color re');
    p.key('Enter');
    await flush();
    expect(p.container.textContent).toContain('Color selection');
    // Takes "red" (the highlighted completion); nothing else is empty (the representations default), so it runs.
    p.key('Enter');
    await flush();
    expect(p.results).toEqual(['Colored the selection red on all']);
  });

  it('shows why a command cannot run', () => {
    const p = mount();
    p.type('measure');
    expect(p.rows()[0]).toContain('Select two atoms first');
  });

  it('selects with a PyMOL expression typed after its name, and explains a broken one', async () => {
    const p = mount();
    p.type('pymol resn HEM and chain A');
    p.key('Enter');
    await flush();
    expect(p.container.textContent).toContain('Select with PyMOL');
    p.key('Enter');
    await flush();
    expect(p.results).toEqual(['Added to the selection: PyMOL "resn HEM and chain A"']);

    const q = mount();
    q.type('pymol (resn HEM');
    q.key('Enter');
    await flush();
    q.key('Enter');
    await flush();
    expect(q.container.querySelector('[role=alert]')?.textContent).toContain('is not closed');
    expect(q.results).toEqual([]);
  });

  it('fills an example expression for editing instead of running it', async () => {
    const p = mount();
    p.type('select with pymol');
    p.key('Enter');
    await flush();
    p.key('ArrowDown'); // the first example: resn HEM
    p.key('Enter');
    await flush();
    expect(p.results).toEqual([]);
    expect(p.input().value).toBe('resn HEM');
    p.type('resn HEM and chain A');
    p.key('Enter');
    await flush();
    expect(p.results).toEqual(['Added to the selection: PyMOL "resn HEM and chain A"']);
  });
});
