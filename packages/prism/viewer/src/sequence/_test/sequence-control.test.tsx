/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, ShellModel } from '@molstar/prism-ui';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { SequenceSource } from '../../sources/sequence.js';
import { SequenceControl } from '../sequence-control.js';
import { sequenceShellUI } from '../shell-ui.js';
import { headlessViewers, hostStructure } from './fixtures.js';

// The band in the page: which section it mounts (a loaded SeqViewSpec's own, sequence.md §6) and how long its tab stays
// (sequence.md §8). jsdom has no 2D canvas: the surface is a stand-in that records the section it was given.

const surfaces = vi.hoisted(() => [] as { section: string; disposed: boolean }[]);
vi.mock('@molstar/seq-viewer/canvas', () => ({
  mountSection: (_ctx: unknown, section: string) => {
    const s = { section, disposed: false };
    surfaces.push(s);
    return { dispose: () => void (s.disposed = true) };
  },
}));
// The control mounts a section only where Seq*'s canvas can draw (it draws through OffscreenCanvas).
vi.stubGlobal('OffscreenCanvas', class {});

afterEach(() => {
  surfaces.length = 0;
  document.body.innerHTML = '';
});

async function until(check: () => boolean, ms = 5000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
}

/** A SeqViewSpec over one structure whose section has a ref of its own (not the built document's `main`). */
function svsOver() {
  const b = createSeqViewSpec({ title: 'Own sections' });
  b.data.complex({ ref: 'own', provider: 'molstar:structure', source: { host_file: '1tup.cif', format: 'mmcif' } });
  const section = b.view({ ref: 'own-view', label: 'Own' }).section({ ref: 'own-main', style: 'compact' });
  section.track({ ref: 'own-letters', label: 'Residues' }).letters({});
  return b.getDocument();
}

const live = () => surfaces.filter((s) => !s.disposed).map((s) => s.section);

describe('the sequence control', () => {
  it("mounts the section the model shows: a loaded SeqViewSpec's own, and the built one's again when it ends", async () => {
    const h = headlessViewers();
    h.structure.set({ ...(await hostStructure()), ref: 'cbs' });
    await until(() => h.viewer.structure.value !== undefined);
    expect(h.viewer.structure.value?.section).toBe('main');
    const host = document.createElement('div');
    document.body.append(host);
    const unmount = domView(() => <SequenceControl viewer={h.viewer} />).mount(host, undefined, undefined as never);
    expect(live()).toEqual(['main']);

    await h.viewer.setOverride({ document: svsOver(), structure: 'cbs' });
    await until(() => h.viewer.structure.value?.section === 'own-main');
    expect(live()).toEqual(['own-main']);
    expect(surfaces.map((s) => s.section)).toEqual(['main', 'own-main']);

    await h.viewer.setOverride(undefined);
    await until(() => h.viewer.structure.value?.section === 'main');
    expect(live()).toEqual(['main']);

    unmount();
    expect(live()).toEqual([]);
    h.viewers.dispose();
  });

  it('keeps its tab across a switch of structure, and drops it once the host shows none', async () => {
    const h = headlessViewers();
    const shell = new ShellModel();
    const registered = vi.spyOn(shell, 'registerTabs');
    const drop = sequenceShellUI({
      viewers: h.viewers,
      source: { structure: h.structure.asReadonly, setBand: () => () => {} } as unknown as SequenceSource,
      shell,
      person: {} as never,
      registerCommands: () => () => {},
      own: (_, parts) => [...parts],
      openCommand: () => {},
      onHover: () => {},
    });
    const tabs: boolean[] = [];
    const sub = shell.dock('top').state.subscribe((s) => tabs.push(s.tabs.includes('sequence')));
    const has = () => shell.dock('top').state.value.tabs.includes('sequence');

    h.structure.set(await hostStructure('1tup'));
    await until(() => h.viewer.structure.value?.title === '1TUP');
    expect(has()).toBe(true);
    tabs.length = 0;
    h.structure.set(await hostStructure('2ocj'));
    expect(h.viewer.structure.value).toBeUndefined();
    await until(() => h.viewer.structure.value?.title === '2OCJ');
    expect(tabs.every((t) => t)).toBe(true);
    expect(registered).toHaveBeenCalledTimes(1);

    h.structure.set(undefined);
    expect(has()).toBe(false);
    sub.unsubscribe();
    drop();
    h.viewers.dispose();
  });
});
