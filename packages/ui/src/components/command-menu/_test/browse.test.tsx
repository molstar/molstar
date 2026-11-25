/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import type { CommandArgs, CommandDefinition } from '../../../commands/types.js';
import { manager, mount, tick } from './fixture.js';

// Browsing and searching (commands §5.1): Recent, groups and categories; a group's shared argument; the highlighted
// row previewed and counted; the keymap's keys; search marks; a blocked row says why; re-ranking.

afterEach(() => {
  document.body.innerHTML = '';
});

function selectGroup() {
  const ran: CommandArgs[] = [];
  const preview = { set: vi.fn(), clear: vi.fn() };
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
  const m = manager([], { preview, browseDelay: 0 });
  m.registerGroups([
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
  m.register([
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
  return { m, ran, preview, t: mount({ manager: m, group: 'select' }) };
}

describe('browsing a group', () => {
  it('shows the shared mode, and runs the commands with it', async () => {
    const { t, ran } = selectGroup();
    const modes = [...t.container.querySelectorAll<HTMLButtonElement>('[aria-label="Select: Mode"] button')];
    expect(modes.map((b) => [b.textContent, b.getAttribute('aria-pressed')])).toEqual([
      ['Replace', 'true'],
      ['Add', 'false'],
    ]);
    modes[1].click();
    await tick();
    expect(modes[1].getAttribute('aria-pressed')).toBe('true');
    t.key('Enter');
    await tick();
    expect(ran).toEqual([{ mode: 'add' }]);
    t.dispose();
  });

  it('previews and counts the highlighted row only, and clears the preview when the menu closes', async () => {
    const { t, m, preview } = selectGroup();
    await tick();
    expect(preview.set).toHaveBeenLastCalledWith('light Ligands');
    expect(t.rows()[0]).toContain('3 residues');
    t.key('ArrowDown');
    await tick();
    expect(preview.set).toHaveBeenLastCalledWith('light Water');
    expect(t.rows()[0]).not.toContain('3 residues');
    expect(t.rows()[1]).toContain('186 residues');
    t.dispose();
    expect(preview.clear).toHaveBeenCalled();
    expect(m.browsing.value).toBeUndefined();
  });

  it('shows the keymap’s key on a row, as caps', async () => {
    const { t, m } = selectGroup();
    m.keymap.remap('other', { key: 'j', alt: true });
    t.type('other');
    await tick();
    const keys = t.container.querySelector('[role=option] kbd');
    expect(keys?.textContent).toMatch(/J$/);
    expect(keys?.querySelectorAll('kbd').length).toBe(2);
    t.dispose();
  });

  it('a crumb names the group; Backspace on an empty input goes up', async () => {
    const { t } = selectGroup();
    expect(t.input().placeholder).toBe('Search in Select');
    t.key('Backspace');
    expect(t.input().placeholder).toBe('Search commands, or type ? for help');
    expect(t.container.textContent).toContain('Browse');
    t.dispose();
  });
});

describe('browsing from the top', () => {
  it('lists Recent, the groups, then the ungrouped commands by category', async () => {
    const m = manager([
      { id: 'a', title: 'Alpha', category: 'View', run: () => ({}) },
      { id: 'b', title: 'Beta', category: 'Look', run: () => ({}) },
      { id: 'c', title: 'Gamma', category: 'Look', group: 'g', run: () => ({}) },
    ]);
    m.registerGroups([{ id: 'g', title: 'Group' }]);
    const r = m.open('b');
    if (r.kind === 'ran') await r.result;
    const t = mount({ manager: m });
    const captions = [...t.container.querySelectorAll('[role=presentation]')].map((c) => c.textContent);
    expect(captions).toEqual(['Recent', 'Browse', 'View']);
    expect(t.rows().map((r) => r.split(/(?=[A-Z])/)[0])).toEqual(['Beta', 'Group', 'Alpha']);
    // → on a group opens it.
    t.key('ArrowDown');
    t.key('ArrowRight');
    expect(t.rows().map((r) => r.slice(0, 5))).toEqual(['Gamma']);
    t.dispose();
  });

  it('marks the matched letters, and a blocked row says why when picked', async () => {
    const m = manager([
      { id: 'a', title: 'Reset the camera', category: 'View', run: () => ({}) },
      { id: 'b', title: 'Remove all', category: 'Edit', when: () => 'Nothing is loaded.', run: () => ({}) },
    ]);
    const t = mount({ manager: m });
    t.type('rem');
    const row = t.options()[0];
    expect(row.getAttribute('aria-disabled')).toBe('true');
    expect([...row.querySelectorAll('mark')].map((x) => x.textContent).join('')).toBe('Rem');
    expect(row.textContent).toContain('Nothing is loaded.');
    t.key('Enter');
    await tick();
    expect(t.container.querySelector('[role=alert]')?.textContent).toBe('Nothing is loaded.');
    expect(t.closed()).toBe(0);
    t.dispose();
  });

  it('shows the middle rows in the order the query ranks them', async () => {
    const titles = ['Color', 'Collapse', 'Outline color', 'Cool', 'Clear color', 'Recolor'];
    const m = manager(titles.map((x) => ({ id: x, title: x, category: 'X', run: () => ({}) })));
    const t = mount({ manager: m });
    const order = (q: string) => m.query(q).matches.map((x) => x.command.title);
    for (const q of ['c', 'co', 'col', 'colo', 'color']) {
      t.type(q);
      await tick();
      const want = order(q);
      expect(t.rows().map((s) => want.find((w) => s.startsWith(w)) ?? s)).toEqual(want);
    }
    t.dispose();
  });
});
