/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { createStorageModel, formatAge, formatBytes, type StorageEntry, type StorageSource } from '../model.js';

// The storage model (journal.md §7): sources registered, listed while watched, Forget and Forget everything.

function fakeSource(id: string, entries: StorageEntry[], o: { forgetAll?: boolean } = {}) {
  const kept = [...entries];
  const source: StorageSource & { kept: StorageEntry[] } = {
    id,
    title: id.toUpperCase(),
    where: 'memory',
    kept,
    list: vi.fn(async () => [...kept]),
    forget: vi.fn(async (entryId: string) => {
      const i = kept.findIndex((e) => e.id === entryId);
      if (i >= 0) kept.splice(i, 1);
    }),
    ...(o.forgetAll
      ? {
          forgetAll: vi.fn(async () => {
            kept.length = 0;
          }),
        }
      : {}),
  };
  return source;
}

describe('the storage model', () => {
  it('registers sources in order, and lists them on refresh with their sizes added up', async () => {
    const model = createStorageModel({ origin: 'https://example.org' });
    expect(model.origin).toBe('https://example.org');
    const a = fakeSource('a', [
      { id: 'x', label: 'X', bytes: 10 },
      { id: 'y', label: 'Y', bytes: 30 },
    ]);
    const b = fakeSource('b', [{ id: 'z', label: 'Z' }]);
    model.register(a);
    model.register(b);
    expect(model.sources.value.map((s) => [s.id, s.listed])).toEqual([
      ['a', false],
      ['b', false],
    ]);
    // Nothing is read until asked (the view is not shown).
    expect(a.list).not.toHaveBeenCalled();
    await model.refresh();
    const [sa, sb] = model.sources.value;
    expect(sa).toMatchObject({ title: 'A', where: 'memory', listed: true, bytes: 40 });
    expect(sa.entries.map((e) => e.id)).toEqual(['x', 'y']);
    // No size known: none said.
    expect(sb.bytes).toBeUndefined();
  });

  it('refuses a second source with the same id, and removes one when its registration is undone', async () => {
    const model = createStorageModel();
    const off = model.register(fakeSource('a', []));
    expect(() => model.register(fakeSource('a', []))).toThrow(/registered already/);
    off();
    expect(model.sources.value).toEqual([]);
    model.register(fakeSource('a', []));
    expect(model.sources.value).toHaveLength(1);
  });

  it('forgets one entry and lists its source again', async () => {
    const model = createStorageModel();
    const a = fakeSource('a', [
      { id: 'x', label: 'X' },
      { id: 'y', label: 'Y' },
    ]);
    model.register(a);
    await model.refresh();
    await model.forget('a', 'x');
    expect(a.forget).toHaveBeenCalledWith('x');
    expect(model.sources.value[0].entries.map((e) => e.id)).toEqual(['y']);
    // An unknown source: nothing happens.
    await model.forget('nope', 'y');
  });

  it('forgets everything: each source by forgetAll where it has one, else entry by entry', async () => {
    const model = createStorageModel();
    const a = fakeSource('a', [{ id: 'x', label: 'X' }], { forgetAll: true });
    const b = fakeSource('b', [
      { id: 'y', label: 'Y' },
      { id: 'z', label: 'Z' },
    ]);
    model.register(a);
    model.register(b);
    await model.forgetEverything();
    expect(a.forgetAll).toHaveBeenCalledTimes(1);
    expect(a.forget).not.toHaveBeenCalled();
    expect(b.forget).toHaveBeenCalledTimes(2);
    expect(model.sources.value.every((s) => s.listed && s.entries.length === 0)).toBe(true);
  });

  it('says when a source cannot be listed or forgotten, and goes on with the others', async () => {
    const model = createStorageModel();
    const bad: StorageSource = {
      id: 'bad',
      title: 'Bad',
      list: async () => {
        throw new Error('blocked');
      },
      forget: async () => {
        throw new Error('refused');
      },
    };
    model.register(bad);
    model.register(fakeSource('a', [{ id: 'x', label: 'X' }]));
    await model.refresh();
    expect(model.sources.value[0].error).toMatch(/blocked/);
    expect(model.sources.value[1].entries).toHaveLength(1);
    await model.forget('bad', 'x');
    expect(model.sources.value[0].error).toMatch(/refused/);
  });

  it('lists when watched, and a source registered while watched at once; not after', async () => {
    const model = createStorageModel();
    const a = fakeSource('a', [{ id: 'x', label: 'X' }]);
    model.register(a);
    const stop = model.watch();
    await vi.waitFor(() => expect(model.sources.value[0].listed).toBe(true));
    const b = fakeSource('b', []);
    model.register(b);
    await vi.waitFor(() => expect(b.list).toHaveBeenCalledTimes(1));
    stop();
    stop();
    const c = fakeSource('c', []);
    model.register(c);
    await Promise.resolve();
    expect(c.list).not.toHaveBeenCalled();
  });

  it('keeps the latest listing when an earlier one ends later', async () => {
    const model = createStorageModel();
    const releases: ((entries: StorageEntry[]) => void)[] = [];
    model.register({
      id: 'slow',
      title: 'Slow',
      list: () => new Promise<StorageEntry[]>((r) => releases.push(r)),
      forget: async () => {},
    });
    const first = model.refresh();
    const second = model.refresh();
    releases[1]([{ id: 'new', label: 'New' }]);
    releases[0]([{ id: 'old', label: 'Old' }]);
    await Promise.all([first, second]);
    expect(model.sources.value[0].entries.map((e) => e.id)).toEqual(['new']);
  });

  it('says sizes and ages in words', () => {
    expect(formatBytes(312)).toBe('312 B');
    expect(formatBytes(4300)).toBe('4.2 KB');
    expect(formatBytes(40_000)).toBe('39 KB');
    expect(formatBytes(1_400_000)).toBe('1.3 MB');
    expect(formatAge(10_000)).toBe('just now');
    expect(formatAge(5 * 60_000)).toBe('5 minutes ago');
    expect(formatAge(60 * 60_000)).toBe('1 hour ago');
    expect(formatAge(3 * 24 * 60 * 60_000)).toBe('3 days ago');
  });
});
