/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterAll, describe, expect, it } from 'vitest';
import 'fake-indexeddb/auto';
import { createMemoryRows, type ThreadRows, type ThreadTable, RowConflict } from '@molstar/ai';
import { createIndexedDbRows } from '../indexeddb-rows.js';

// The row store contract, directly, on each backend: what the host relies on (analysis/viewer-mode.md §6).
const closers: (() => void)[] = [];
afterAll(() => {
  for (const close of closers) close();
});
let n = 0;
const closing = <S extends { close(): void }>(store: S): S => {
  closers.push(() => store.close());
  return store;
};

interface Backend {
  readonly name: string;
  /** Two handles to one fresh storage. */
  open(): [ThreadRows, ThreadRows];
  readonly watches: boolean;
}
const BACKENDS: Backend[] = [
  {
    name: 'memory',
    open: () => {
      const rows = createMemoryRows<ThreadTable>();
      return [rows, rows];
    },
    watches: false,
  },
  {
    name: 'indexeddb',
    open: () => {
      const name = `rows-${++n}`;
      return [closing(createIndexedDbRows<ThreadTable>(name)), closing(createIndexedDbRows<ThreadTable>(name))];
    },
    watches: true,
  },
];

const put = (key: string, json: string) => ({ table: 'entry' as const, key, row: { json } });
const wait = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));

for (const backend of BACKENDS)
  describe(`thread rows: ${backend.name}`, () => {
    it('has no rows and generation 0 for a scope never written', async () => {
      const [rows] = backend.open();
      expect(await rows.load('s')).toBeUndefined();
      expect(await rows.generation('s')).toBe(0);
    });

    it('applies a write set and moves the generation by one, deleting what is written without a row', async () => {
      const [rows] = backend.open();
      await rows.apply('s', 0, [put('a', '1'), put('b', '2'), { table: 'thread', key: 'thread', row: { json: '{}' } }]);
      await rows.apply('s', 1, [put('a', '3'), { table: 'entry', key: 'b', row: undefined }]);
      const loaded = await rows.load('s');
      expect(loaded?.generation).toBe(2);
      expect(await rows.generation('s')).toBe(2);
      expect([...(loaded?.rows ?? [])].sort((a, b) => (a.table + a.key).localeCompare(b.table + b.key))).toEqual([
        put('a', '3'),
        { table: 'thread', key: 'thread', row: { json: '{}' } },
      ]);
    });

    it('keeps scopes apart, and drops one with its generation', async () => {
      const [rows] = backend.open();
      await rows.apply('s', 0, [put('a', '1')]);
      await rows.apply('t', 0, [put('a', '2')]);
      await rows.drop('s');
      expect(await rows.load('s')).toBeUndefined();
      expect(await rows.generation('s')).toBe(0);
      expect((await rows.load('t'))?.rows).toEqual([put('a', '2')]);
    });

    it('refuses a write against a stale generation, whole, whichever handle wrote first', async () => {
      const [first, second] = backend.open();
      await first.apply('s', 0, [put('a', '1')]);
      await expect(second.apply('s', 0, [put('b', '2'), put('c', '3')])).rejects.toBeInstanceOf(RowConflict);
      expect((await second.load('s'))?.rows).toEqual([put('a', '1')]);
      expect(await second.generation('s')).toBe(1);
      await second.apply('s', 1, [put('b', '2')]);
      expect((await first.load('s'))?.generation).toBe(2);
    });

    it('lets exactly one of simultaneous writes at the same generation through', async () => {
      const [first, second] = backend.open();
      const results = await Promise.allSettled([
        first.apply('s', 0, [put('a', 'first')]),
        second.apply('s', 0, [put('a', 'second')]),
      ]);
      expect(results.map((r) => r.status).sort()).toEqual(['fulfilled', 'rejected']);
      expect(await first.generation('s')).toBe(1);
    });

    if (backend.watches)
      it('tells the other handle which scope was written or dropped, and not the one that wrote', async () => {
        const [first, second] = backend.open();
        const heard: string[] = [];
        const own: string[] = [];
        const stop = second.watch?.((scope) => heard.push(scope));
        const stopOwn = first.watch?.((scope) => own.push(scope));
        await first.apply('s', 0, [put('a', '1')]);
        await first.drop('t');
        await wait(50);
        stop?.();
        stopOwn?.();
        await first.apply('s', 1, [put('a', '2')]);
        await wait(50);
        expect({ heard, own }).toEqual({ heard: ['s', 't'], own: [] });
      });
  });

it('indexeddb refuses a database of the same name and version made by something else', async () => {
  const other = indexedDB.open('not-threads', 1);
  other.onupgradeneeded = () => other.result.createObjectStore('agent');
  await new Promise((resolve) => (other.onsuccess = resolve));
  other.result.close();
  const rows = closing(createIndexedDbRows<ThreadTable>('not-threads'));
  await expect(rows.generation('s')).rejects.toThrow('IndexedDB database "not-threads" holds something else');
});
