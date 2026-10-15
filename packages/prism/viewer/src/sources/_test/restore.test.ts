/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { memoryStore } from '../../persist/stores.js';

// The light restore (journal.md §6) on headless viewers sharing one store, as a tab's sessionStorage is across a reload.
// The address is served by a stubbed fetch: no network.

const cif = fs.readFileSync(
  path.resolve(import.meta.dirname, '../../../../../../data/examples/1cbs_updated.cif'),
  'utf8',
);
const Url = 'https://files.rcsb.org/download/1cbs.cif';

afterEach(() => {
  vi.unstubAllGlobals();
});

const structures = (v: {
  plugin: { managers: { structure: { hierarchy: { current: { structures: unknown[] } } } } };
}) => v.plugin.managers.structure.hierarchy.current.structures.length;

describe('the light restore', () => {
  it('keeps a scene loaded from an address and brings it back after a reload, with an empty history', async () => {
    vi.stubGlobal('fetch', async () => new Response(cif));
    const store = memoryStore({ version: 'test', budget: 256 * 1024 });
    const first = await headlessViewer({ spec: { persist: { key: 'test', restore: store } } });
    expectOk(await first.api.load.url!({ url: Url, format: 'mmcif' }));
    expectOk(await first.api.styles.apply({ style: 'surface' }));
    first.restore!.flush();
    expect(first.restore!.available).toBe(true);

    const second = await headlessViewer({ spec: { persist: { key: 'test', restore: store } } });
    expect(structures(second)).toBe(0);
    expect(await second.restore!.run()).toBe(true);
    expect(structures(second)).toBe(1);
    expect(second.managers.history.steps.value).toEqual([]);
  });

  it('does not keep a scene with a local file, and says so once', async () => {
    const store = memoryStore({ version: 'test', budget: 256 * 1024 });
    const v = await headlessViewer({ spec: { persist: { key: 'test', restore: store } } });
    const said: string[] = [];
    v.managers.notifications.posted.subscribe((m) => said.push(m.text));
    expectOk(await v.api.load.files({ files: [cbs()] }));
    v.restore!.flush();
    v.restore!.flush();
    expect(v.restore!.available).toBe(false);
    expect(said.filter((t) => /not kept across a reload/.test(t)).length).toBeGreaterThan(0);
  });

  it('drops what another version wrote, and anything over the budget', () => {
    const a = memoryStore({ version: '1', budget: 100 });
    a.set('x', { small: true });
    expect(a.get('x')).toEqual({ small: true });
    a.set('x', { big: 'y'.repeat(200) });
    expect(a.get('x')).toBeUndefined();
  });
});
