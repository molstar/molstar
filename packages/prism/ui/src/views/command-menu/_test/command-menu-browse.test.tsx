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
import { CommandManager } from '../../../commands/manager.js';
import type { CommandArgs, CommandDefinition } from '../../../commands/types.js';
import { render } from '../../../dom/render.js';
import { flushFrame } from '../../../dom/scope.js';
import { observableState } from '../../../reactive/observable-state.js';
import { CommandMenu } from '../command-menu.js';

// Browsing a group with a shared argument (commands.md §5): the mode in the header, previews and counts for the
// highlighted row, the keymap's keys on the rows.

const tick = () => new Promise((r) => setTimeout(r, 0));

afterEach(() => {
  document.body.innerHTML = '';
});

function setup() {
  const ran: CommandArgs[] = [];
  const preview = { set: vi.fn(), clear: vi.fn() };
  const manager = new CommandManager<object>({
    context: observableState<object>({}).asReadonly,
    formats: { get: () => undefined },
    preview,
    browseDelay: 0,
  });
  const query = (id: string, n: number): CommandDefinition<object> => ({
    id,
    title: id,
    category: 'Select',
    group: 'select',
    preview: () => `light ${id}`,
    count: async () => `${n} residues`,
    run: (args) => {
      ran.push(args);
      return { message: id };
    },
  });
  manager.registerGroups([
    {
      id: 'select',
      title: 'Select',
      shared: {
        name: 'mode',
        title: 'Mode',
        options: [
          ['set', 'Replace'],
          ['add', 'Add'],
        ],
        default: 'set',
      },
    },
  ]);
  manager.register([
    query('Ligands', 3),
    query('Water', 186),
    {
      id: 'other',
      title: 'Other',
      category: 'Misc',
      group: 'select',
      shortcut: { key: 'o', alt: true },
      run: () => ({}),
    },
  ]);
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <CommandMenu manager={manager} group="select" onClose={() => {}} />, container);
  const input = () => container.querySelector('input') as HTMLInputElement;
  const key = (k: string) => input().dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  const rows = () => [...container.querySelectorAll('[role=option]')].map((r) => r.textContent ?? '');
  return { manager, preview, ran, container, key, rows, dispose };
}

describe('browsing a group', () => {
  it('shows the shared mode, and runs the commands with it', async () => {
    const p = setup();
    const modes = [...p.container.querySelectorAll<HTMLButtonElement>('[aria-label="Select: Mode"] button')];
    expect(modes.map((b) => [b.textContent, b.getAttribute('aria-pressed')])).toEqual([
      ['Replace', 'true'],
      ['Add', 'false'],
    ]);
    modes[1]!.click();
    await tick();
    expect(modes[1]!.getAttribute('aria-pressed')).toBe('true');
    p.key('Enter');
    await tick();
    expect(p.ran).toEqual([{ mode: 'add' }]);
    p.dispose();
  });

  it('previews and counts the highlighted row only, and clears the preview when the menu closes', async () => {
    const p = setup();
    await tick();
    expect(p.preview.set).toHaveBeenLastCalledWith('light Ligands');
    expect(p.rows()[0]).toContain('3 residues');
    p.key('ArrowDown');
    await tick();
    expect(p.preview.set).toHaveBeenLastCalledWith('light Water');
    expect(p.rows()[0]).not.toContain('3 residues');
    expect(p.rows()[1]).toContain('186 residues');
    p.dispose();
    expect(p.preview.clear).toHaveBeenCalled();
    expect(p.manager.browsing.value).toBeUndefined();
  });

  it('shows the keymap’s key on a row', async () => {
    const p = setup();
    p.manager.keymap.remap('other', { key: 'j', alt: true });
    input(p).value = 'other';
    input(p).dispatchEvent(new Event('input'));
    // The frame the input asks for: the list is filtered once a frame (rendering.md §2.3).
    flushFrame();
    await tick();
    expect(p.container.querySelector('[role=option] kbd')?.textContent).toMatch(/J$/);
    p.dispose();
  });
});

const input = (p: { container: HTMLElement }) => p.container.querySelector('input') as HTMLInputElement;
