/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';

// The Select group's commands are made when wanted (commands/select.ts, commands.md §4 `registerLazy`): loading N
// structures one by one lists the queries over the loaded structures not at all (it was once per load, over every
// structure: N²); the first read lists them once, over all of them.

const lists = vi.hoisted(() => ({ n: 0 }));
vi.mock('../../sources/select-queries.js', async (importOriginal) => {
  const actual = await importOriginal<typeof import('../../sources/select-queries.js')>();
  return {
    ...actual,
    selectQuerySource: (plugin: Parameters<typeof actual.selectQuerySource>[0]) => {
      const source = actual.selectQuerySource(plugin);
      return {
        ...source,
        list: () => {
          lists.n++;
          return source.list();
        },
      };
    },
  };
});

describe('the Select group, made when wanted', () => {
  it('lists nothing while structures load; the first read lists once; the menu follows changes', async () => {
    const v = await headlessViewer();
    try {
      const commands = v.managers.commands;
      commands.get('select.set');
      lists.n = 0;
      for (let i = 0; i < 4; i++) expectOk(await v.api.load.files({ files: [cbs()] }));
      expect(lists.n).toBe(0);
      expect(commands.get('select.chain.A')?.title).toBe('Chain A');
      commands.query('chain');
      expect(lists.n).toBe(1);

      // Open (the menu holds `want`), a change is listed in a microtask, not by the menu's reads.
      const release = commands.want();
      expectOk(await v.api.load.files({ files: [cbs()] }));
      await Promise.resolve();
      expect(lists.n).toBe(2);
      release();
    } finally {
      await v.dispose();
    }
  });

  it('a key bound to a structure query command runs it after a load, the menu never opened', async () => {
    const v = await headlessViewer();
    try {
      const commands = v.managers.commands;
      commands.keymap.setApp({ 'select.chain.A': { key: 'j' } });
      expectOk(await v.api.load.files({ files: [cbs()] }));
      lists.n = 0;
      expect(commands.handleKey({ key: 'j' })).toBe(true);
      expect(lists.n).toBe(1);
    } finally {
      await v.dispose();
    }
  });
});
