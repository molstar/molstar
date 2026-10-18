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
import { mockCommands, mockPrefixes } from '../testing/commands.js';
import { mockSnapshots } from '@molstar/prism-ui/testing/index';
import { mockScene } from '../testing/scene.js';
import { render } from '@molstar/prism-ui/dom/render';
import { flushFrame } from '@molstar/prism-ui/dom/scope';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { CommandMenu } from '@molstar/prism-ui';

const flush = () => new Promise((r) => setTimeout(r, 0));

function mount(options: { prefixes?: boolean } = {}) {
  const shell = new ShellModel();
  const manager = mockCommands(shell);
  const prefixes = options.prefixes
    ? mockPrefixes({ manager, shell, scene: mockScene(), snapshots: mockSnapshots(() => {}) })
    : undefined;
  const closed = vi.fn();
  const results: string[] = [];
  manager.outcomes.subscribe((o) => {
    if (o.kind === 'done') results.push(o.result.message ?? '');
  });
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <CommandMenu manager={manager} prefixes={prefixes} onClose={closed} />, container);
  const input = () => container.querySelector('input') as HTMLInputElement;
  /** Typed text, without the frame that filters it (`type` runs it). */
  const typeOnly = (text: string) => {
    input().value = text;
    input().dispatchEvent(new Event('input'));
  };
  const type = (text: string) => {
    typeOnly(text);
    // The frame the input asks for: the list is filtered once a frame (rendering.md §2.3).
    flushFrame();
  };
  const key = (k: string) => input().dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  const rows = () => [...container.querySelectorAll('[role=option]')].map((r) => r.textContent ?? '');
  return { manager, closed, results, container, input, type, typeOnly, key, rows, dispose };
}

afterEach(() => {
  document.body.innerHTML = '';
});

describe('CommandMenu', () => {
  it('matches once per keystroke, and not at all once closed', () => {
    const p = mount();
    const query = vi.spyOn(p.manager, 'query');
    // The list and the highlight's bounds read one computation of the rows.
    p.type('reset zo');
    expect(query).toHaveBeenCalledTimes(1);
    p.dispose();
    query.mockClear();
    p.manager.register([{ id: 'test.late', title: 'Late', category: 'Test', run: () => {} } as never]);
    expect(query).not.toHaveBeenCalled();
  });

  it('filters once for keys typed within a frame, for the latest text only', () => {
    const p = mount();
    const query = vi.spyOn(p.manager, 'query');
    // Ten characters faster than frames (a fast typist, a held key): each one shown in the input, none filtered yet.
    const text = 'reset zoom';
    for (let i = 1; i <= text.length; i++) p.typeOnly(text.slice(0, i));
    expect(p.input().value).toBe(text);
    expect(query).not.toHaveBeenCalled();
    flushFrame();
    expect(query.mock.calls.map((c) => c[0])).toEqual([text]);
    expect(p.rows()[0]).toContain('Reset zoom');
    p.dispose();
  });

  it('a key acting on the rows filters the latest text first', async () => {
    const p = mount();
    p.typeOnly('reset zo');
    // Enter before the frame: the rows it picks from are the latest text's.
    p.key('Enter');
    await flush();
    expect(p.results).toEqual(['Zoom reset']);
    p.dispose();
  });

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

  it('browses groups: opens them, shows the path, goes back up', () => {
    const p = mount();
    const row = (text: string) =>
      [...p.container.querySelectorAll('[role=option]')].findIndex((r) => r.textContent?.startsWith(text));
    // Down to "Select", open it; then "Expression".
    for (let i = 0; i < row('Select'); i++) p.key('ArrowDown');
    p.key('Enter');
    expect(p.container.querySelector('[class*="crumb"]')?.textContent).toBe('Select');
    for (let i = 0; i < row('Expression'); i++) p.key('ArrowDown');
    p.key('ArrowRight');
    expect(p.rows().map((r) => r.split('…')[0])).toEqual(['Select with PyMOL', 'Select with VMD', 'Select with Jmol']);
    p.key('Backspace');
    expect(p.container.querySelector('[class*="crumb"]')?.textContent).toBe('Select');
    // A search shows a grouped command's path.
    p.key('Backspace');
    p.type('pymol');
    expect(p.rows()[0]).toContain('Select › Expression');
  });

  it('turns into another list with a prefix: ":" opens a select form, "?" shows help', async () => {
    const p = mount({ prefixes: true });
    p.type(':pymol resn HEM');
    expect(p.rows()[0]).toContain('Select with PyMOL: resn HEM');
    p.key('Enter');
    await flush();
    expect(p.manager.active.value?.command.id).toBe('select.pymol');
    expect(p.closed).not.toHaveBeenCalled();

    const q = mount({ prefixes: true });
    q.type('?');
    expect(q.rows().join('|')).toContain('@  Go to');
    expect(q.rows().join('|')).toContain('Toggle illumination');
  });

  it('finds commands by their aliases (align → Superpose, screenshot → Export image)', () => {
    const shell = new ShellModel();
    const manager = mockCommands(shell, () => 2);
    const first = (text: string) => manager.query(text).matches[0]?.command.id;
    expect(first('align')).toBe('superpose.open');
    expect(first('screenshot')).toBe('export.image');
    expect(first('deselect')).toBe('selection.clear');
  });
});
