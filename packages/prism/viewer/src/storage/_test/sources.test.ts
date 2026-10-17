/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import 'fake-indexeddb/auto';
import { createIndexedDbRows } from '../../journal/indexeddb-rows.js';
import { createMemoryRows, type ThreadTable } from '@molstar/ai';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { fakeApi } from '../../journal/_test/runs.js';
import { createRecordIndex, createSessionRecord, keptSize, RecordsDatabase } from '../../journal/record.js';
import { memoryStore, webStore } from '../../persist/stores.js';
import { keptRecordsSource, storeSource, viewerStorageSources } from '../sources.js';

// The viewer's own storage sources (journal.md §7): kept records (memory rows, and IndexedDB read again as a new page
// would), the preferences and the light restore.

const DayMs = 24 * 60 * 60_000;

class MemoryStorage {
  readonly map = new Map<string, string>();
  get length() {
    return this.map.size;
  }
  key(i: number) {
    return [...this.map.keys()][i] ?? null;
  }
  getItem(k: string) {
    return this.map.get(k) ?? null;
  }
  setItem(k: string, v: string) {
    this.map.set(k, v);
  }
  removeItem(k: string) {
    this.map.delete(k);
  }
}

afterEach(() => {
  vi.unstubAllGlobals();
});

/** A record kept in these rows at `at`, with a change or two. */
async function keptRecord(rows: ReturnType<typeof createMemoryRows<ThreadTable>>, title: string, at: number) {
  const { api, run } = fakeApi();
  const record = await createSessionRecord({ api, clock: () => at });
  run();
  await record.keep(rows, await createRecordIndex(rows), title);
  run();
  await record.flush();
  return record;
}

describe('kept records', () => {
  it('lists each with its title, size and age; marks this session and those older than 90 days', async () => {
    const rows = createMemoryRows<ThreadTable>();
    const t0 = 1_000_000_000_000;
    const old = await keptRecord(rows, '1CBS', t0);
    await old.dispose();
    const current = await keptRecord(rows, '4HHB', t0 + 100 * DayMs);
    const source = keptRecordsSource({ rows, record: current, now: () => t0 + 100 * DayMs + 5 * 60_000 });
    const entries = await source.list();
    expect(entries.map((e) => e.label)).toEqual(['4HHB', '1CBS']);
    const [now, then] = entries;
    expect(now.detail).toBe('This session · last written 5 minutes ago');
    expect(now.stale).toBeUndefined();
    expect(now.bytes).toBeGreaterThan(0);
    expect(then.stale).toBe(true);
    expect(then.detail).toMatch(/^Last written 100 days ago · older than 90 days/);
    await current.dispose();
  });

  it('forgets another record: its rows and its index entry', async () => {
    const rows = createMemoryRows<ThreadTable>();
    const old = await keptRecord(rows, '1CBS', Date.now());
    await old.dispose();
    const source = keptRecordsSource({ rows });
    expect(await keptSize(rows, old.id)).toBeGreaterThan(0);
    await source.forget(old.id);
    expect(await source.list()).toEqual([]);
    expect(await keptSize(rows, old.id)).toBeUndefined();
  });

  it("forgets this session's kept record: it goes on in memory with what it holds, and can be kept again", async () => {
    const rows = createMemoryRows<ThreadTable>();
    const record = await keptRecord(rows, '1CBS', Date.now());
    const before = (await record.read(10)).items.map((i) => i.id);
    expect(before).toHaveLength(2);
    const source = keptRecordsSource({ rows, record });
    await source.forget(record.id);
    expect(record.kept).toBe(false);
    expect(await source.list()).toEqual([]);
    expect(await keptSize(rows, record.id)).toBeUndefined();
    // History still reads it, and it goes on writing in memory.
    expect((await record.read(10)).items.map((i) => i.id)).toEqual(before);
    await record.keep(rows, await createRecordIndex(rows), '1CBS again');
    expect(record.kept).toBe(true);
    expect((await source.list()).map((e) => e.label)).toEqual(['1CBS again']);
    await record.dispose();
  });

  it('reads IndexedDB only when the database is there, and lists what a page before kept', async () => {
    const source = keptRecordsSource();
    expect(await source.list()).toEqual([]);
    // Listing made no empty database.
    expect((await indexedDB.databases()).map((d) => d.name)).not.toContain(RecordsDatabase);

    const rows = createIndexedDbRows<ThreadTable>(RecordsDatabase);
    const record = await keptRecord(rows, '1CBS', Date.now());
    await record.dispose();
    rows.close();

    const entries = await source.list();
    expect(entries.map((e) => e.label)).toEqual(['1CBS']);
    expect(entries[0].bytes).toBeGreaterThan(0);
    await source.forget(entries[0].id);
    expect(await source.list()).toEqual([]);
    source.dispose();
  });

  it('lists nothing where there is no IndexedDB', async () => {
    vi.stubGlobal('indexedDB', undefined);
    expect(await keptRecordsSource().list()).toEqual([]);
  });
});

describe('the light stores', () => {
  it('lists a store’s values with labels and sizes; Forget removes one', async () => {
    const store = memoryStore({ version: '1', budget: 1000 });
    store.set('theme', 'dark');
    store.set('other', { a: 1 });
    const source = storeSource({
      id: 'preferences',
      title: 'Preferences',
      where: 'memory',
      store,
      labels: { theme: 'Theme' },
    });
    const entries = await source.list();
    expect(entries.map((e) => e.label).sort()).toEqual(['Theme', 'other']);
    expect(entries.every((e) => (e.bytes ?? 0) > 0)).toBe(true);
    await source.forget('theme');
    expect((await source.list()).map((e) => e.id)).toEqual(['other']);
    expect(store.get('theme')).toBeUndefined();
  });

  it('drops values of another version when listing, and leaves what is not its own', async () => {
    const local = new MemoryStorage();
    vi.stubGlobal('localStorage', local);
    const store = webStore('local', { prefix: 'prism-viewer', version: '2', budget: 1000 });
    store.set('theme', 'dark');
    local.setItem('prism-viewer.layout', JSON.stringify({ v: '1', d: {} }));
    local.setItem('prism-viewer.notours', 'plain text');
    local.setItem('prism-assistant.setup', JSON.stringify({ v: '1', d: {} }));
    const entries = store.entries();
    expect(entries.map((e) => e.key)).toEqual(['theme']);
    expect(entries[0].bytes).toBe((local.getItem('prism-viewer.theme')!.length + 'prism-viewer.theme'.length) * 2);
    expect(local.getItem('prism-viewer.layout')).toBeNull();
    expect(local.getItem('prism-viewer.notours')).toBe('plain text');
    expect(local.getItem('prism-assistant.setup')).not.toBeNull();
  });

  it('Forget on the keyboard shortcuts resets the keymap’s in-memory layer, and nothing else does', async () => {
    const preferences = memoryStore({ version: '1', budget: 1000 });
    preferences.set('keybindings', { a: { key: 'j' } });
    preferences.set('theme', 'dark');
    const reset = vi.fn();
    const prefs = viewerStorageSources({ persist: { key: 'v', preferences }, resetKeybindings: reset }).sources[1];
    await prefs.forget('theme');
    expect(reset).not.toHaveBeenCalled();
    await prefs.forget('keybindings');
    expect(reset).toHaveBeenCalledTimes(1);
    expect(preferences.get('keybindings')).toBeUndefined();
  });

  it("are the viewer's sources when the spec gives listing stores", async () => {
    const local = new MemoryStorage();
    const session = new MemoryStorage();
    vi.stubGlobal('localStorage', local);
    vi.stubGlobal('sessionStorage', session);
    const preferences = webStore('local', { prefix: 'v', version: '1', budget: 1000 });
    const restore = webStore('session', { prefix: 'v', version: '1', budget: 1000 });
    preferences.set('theme', 'light');
    restore.set('restore', { snapshot: {} });
    const own = viewerStorageSources({ persist: { key: 'v', preferences, restore } });
    expect(own.sources.map((s) => s.id)).toEqual(['kept-records', 'preferences', 'light-restore']);
    const [, prefs, light] = own.sources;
    expect((await prefs.list()).map((e) => e.label)).toEqual(['Theme']);
    expect((await light.list()).map((e) => e.label)).toEqual(['This tab’s scene, without its data']);
    await light.forget('restore');
    expect(session.getItem('v.restore')).toBeNull();
    own.dispose();
    // Stores that do not list (a host's own): only the kept records.
    const plain = { get: () => undefined, set: () => {}, remove: () => {} };
    expect(viewerStorageSources({ persist: { key: 'v', preferences: plain } }).sources.map((s) => s.id)).toEqual([
      'kept-records',
    ]);
  });
});
