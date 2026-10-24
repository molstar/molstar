/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CanvasMenuItem, ShellModel } from '@molstar/prism-ui';
import { firstValueFrom } from 'rxjs';
import { describe, expect, it } from 'vitest';
import { sequenceCommands, sequenceMenu } from '../commands.js';
import { headlessViewers } from './fixtures.js';

const tab = {
  id: 'sequence',
  title: 'Sequence',
  icon: 'sequence',
  dock: 'top',
  frame: 'none',
  render: { mount: () => () => {} },
} as const;

describe('the Sequence menu', () => {
  const setup = () => {
    const h = headlessViewers();
    const shell = new ShellModel();
    shell.registerTabs([tab as never]);
    shell.show('sequence');
    const commands = sequenceCommands(h.viewer, shell);
    const run = (id: string) => commands.find((c) => c.id === id)?.run({}, {}, new AbortController().signal);
    const open = (id: string) => void run(id);
    const menu = sequenceMenu(h.viewer, shell, open);
    const marked = async () => {
      const out: string[] = [];
      for (const m of menu) if (m.choice && m.pressed && (await firstValueFrom(m.pressed))) out.push(m.id);
      return out;
    };
    return { h, shell, run, menu, marked };
  };
  const item = (menu: CanvasMenuItem[], id: string) => menu.find((m) => m.id === id) as CanvasMenuItem;

  it('is one single choice of No sequence, Text and Tracks (Tracks the default), each with its icon', async () => {
    const { menu, marked, h } = setup();
    const view = menu.filter((m) => m.group === 'view');
    expect(view.map((m) => [m.id, m.icon, m.choice])).toEqual([
      ['none', 'ban', true],
      ['text', 'align-left', true],
      ['tracks', 'chart-gantt', true],
    ]);
    expect((await marked()).filter((id) => ['none', 'text', 'tracks'].includes(id))).toEqual(['tracks']);
    h.viewers.dispose();
  });

  it('picking one hides the band, or shows it in that style', async () => {
    const { menu, marked, shell, run, h } = setup();
    item(menu, 'text').run();
    expect(h.viewer.style.value).toBe('compact');
    expect(shell.dock('top').active.value?.id).toBe('sequence');
    expect((await marked()).filter((id) => ['none', 'text', 'tracks'].includes(id))).toEqual(['text']);
    item(menu, 'none').run();
    expect(shell.dock('top').active.value).toBeUndefined();
    expect((await marked()).filter((id) => ['none', 'text', 'tracks'].includes(id))).toEqual(['none']);
    // From hidden, a style shows it again.
    item(menu, 'tracks').run();
    expect(shell.dock('top').active.value?.id).toBe('sequence');
    expect(h.viewer.style.value).toBe('tracks');
    // The rail is a command, not a menu item.
    expect(menu.some((m) => m.id.startsWith('rail'))).toBe(false);
    await run('sequence.rail.strip');
    expect(h.viewer.rail.value.navigator).toBe('strip');
    h.viewers.dispose();
  });

  it('holds only No sequence, Text and Tracks, then Edit in structure tree where there is a tree', async () => {
    const { menu, h, shell } = setup();
    expect(menu.map((m) => m.id)).toEqual(['none', 'text', 'tracks']);
    expect(sequenceMenu(h.viewer, shell, () => {}, true).map((m) => m.id)).toEqual([
      'none',
      'text',
      'tracks',
      'reveal',
    ]);
    h.viewers.dispose();
  });
});
