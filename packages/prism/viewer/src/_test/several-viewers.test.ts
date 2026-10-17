/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import 'fake-indexeddb/auto';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { keepSessionRecord } from '../ui/history.js';
import { createPrismViewer } from '../page.js';
import { webStore } from '../persist/stores.js';
import { VmdLanguage } from '../targets/languages/vmd.js';
import type { PrismViewer } from '../viewer.js';
import { recordsDatabase } from '../journal/record.js';

// Several real viewers on one page (README decision 29), over jsdom (no WebGL: their canvases do not start): what each
// keeps apart. Keys, menus and the pick level are covered where they live (prism-ui's `ownsKey` tests, the
// headless two-viewer tests in sources/_test/restore.test.ts and element.test.ts).

// No WebGL in jsdom: the canvas does not start (the shell says so).
vi.mock('../sources/canvas-host.js', () => ({ mountCanvas: async () => false, unmountCanvas() {} }));

class NoResize {
  observe() {}
  unobserve() {}
  disconnect() {}
}

const open: PrismViewer[] = [];
beforeEach(() => {
  vi.stubGlobal('ResizeObserver', NoResize);
  vi.spyOn(console, 'error').mockImplementation(() => {});
});
afterEach(async () => {
  for (const v of open.splice(0)) await v.dispose();
  document.body.innerHTML = '';
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

async function make(key: string, spec: Partial<ReturnType<typeof DefaultPrismViewerSpec>> = {}) {
  const target = document.createElement('div');
  document.body.append(target);
  const viewer = await createPrismViewer({
    target,
    spec: {
      ...DefaultPrismViewerSpec(),
      persist: { key, preferences: webStore('local', { prefix: key, version: '1', budget: 64 * 1024 }) },
      ...spec,
    },
  });
  open.push(viewer);
  return { viewer, target };
}

describe('several viewers on a page', () => {
  it('refuse a second with a persist key in use, and free it when the first goes', async () => {
    const a = await make('one');
    await expect(make('one')).rejects.toThrow(/keeps its state under 'one'/);
    await a.viewer.dispose();
    open.splice(open.indexOf(a.viewer), 1);
    const again = await make('one');
    expect(again.viewer).toBeDefined();
  });

  it('keep their own preferences under their own keys', async () => {
    const a = await make('left');
    const b = await make('right');
    a.viewer.managers.commands.open('view.theme');
    await new Promise((r) => setTimeout(r, 0));
    const keys = Object.keys(localStorage);
    expect(keys.some((k) => k.startsWith('left'))).toBe(true);
    expect(keys.some((k) => k.startsWith('right'))).toBe(false);
    expect(b.viewer.spec.persist?.key).toBe('right');
  });

  it('read their own query languages', async () => {
    const a = await make('qa', { languages: [VmdLanguage] });
    const b = await make('qb', { languages: [] });
    expect(a.viewer.api.describe('select.set')?.params.language).toMatchObject({ options: [['vmd', 'VMD']] });
    expect(b.viewer.api.describe('select.set')?.params.language).toMatchObject({ options: [] });
  });

  it('go fullscreen one at a time, each following only its own element', async () => {
    const a = await make('fa');
    const b = await make('fb');
    const enter = (target: HTMLElement | null) => {
      Object.defineProperty(document, 'fullscreenElement', { value: target, configurable: true });
      (target ?? a.target).dispatchEvent(new Event('fullscreenchange', { bubbles: true }));
    };
    try {
      enter(a.target);
      expect([a.viewer.managers.navigation.fullscreen.value, b.viewer.managers.navigation.fullscreen.value]).toEqual([
        true,
        false,
      ]);
      enter(b.target);
      expect([a.viewer.managers.navigation.fullscreen.value, b.viewer.managers.navigation.fullscreen.value]).toEqual([
        false,
        true,
      ]);
      enter(null);
      expect([a.viewer.managers.navigation.fullscreen.value, b.viewer.managers.navigation.fullscreen.value]).toEqual([
        false,
        false,
      ]);
    } finally {
      delete (document as unknown as Record<string, unknown>).fullscreenElement;
    }
  });

  it('keep their records in databases of their own', async () => {
    expect(recordsDatabase('fa')).not.toBe(recordsDatabase('fb'));
    // A viewer without a key never shares the app's ('prism-viewer'); unkeyed viewers share theirs.
    expect(recordsDatabase()).toBe('prism-default-records');
    expect(recordsDatabase('')).toBe('prism-default-records');
    expect(recordsDatabase('prism-viewer')).toBe('prism-viewer-records');

    const a = await make('ra');
    const b = await make('rb');
    expect(await keepSessionRecord(a.viewer)).toBe(true);
    const names = (await indexedDB.databases()).map((d) => d.name);
    expect(names).toContain('ra-records');
    expect(names).not.toContain('rb-records');
    expect(b.viewer.record.kept).toBe(false);

    expect(await keepSessionRecord(b.viewer)).toBe(true);
    expect((await indexedDB.databases()).map((d) => d.name)).toEqual(
      expect.arrayContaining(['ra-records', 'rb-records']),
    );

    // Each lists only its own kept records (the other's is not here).
    await a.viewer.keptRecords.refresh();
    await b.viewer.keptRecords.refresh();
    expect(a.viewer.keptRecords.list.value.map((r) => r.id)).not.toContain(b.viewer.record.id);
    expect(b.viewer.keptRecords.list.value.map((r) => r.id)).not.toContain(a.viewer.record.id);
  });
});

describe('the light stores of several viewers', () => {
  it('keep apart by prefix', () => {
    const a = webStore('local', { prefix: 'pa', version: '1', budget: 1000 });
    const b = webStore('local', { prefix: 'pb', version: '1', budget: 1000 });
    a.set('k', 1);
    expect(b.get('k')).toBeUndefined();
  });
});
