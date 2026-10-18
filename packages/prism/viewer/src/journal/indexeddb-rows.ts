/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Row, RowConflict, type RowStore, type RowWrite } from '@molstar/ai';

const ROWS = 'rows';
const GENERATIONS = 'generations';
const STORES = [ROWS, GENERATIONS] as const;

interface RowRecord {
  scope: string;
  table: string;
  key: string;
  row: Row;
}

const request = <R>(r: IDBRequest<R>): Promise<R> =>
  new Promise((resolve, reject) => {
    r.onsuccess = () => resolve(r.result);
    r.onerror = () => reject(r.error);
  });

const done = (tx: IDBTransaction): Promise<void> =>
  new Promise((resolve, reject) => {
    tx.oncomplete = () => resolve();
    tx.onerror = () => reject(tx.error);
    tx.onabort = () => reject(tx.error ?? new Error('transaction aborted'));
  });

const scopeRange = (scope: string) => IDBKeyRange.bound([scope], [scope, []]);

function open(name: string): Promise<IDBDatabase> {
  const opening = indexedDB.open(name, 1);
  opening.onupgradeneeded = () => {
    opening.result.createObjectStore(ROWS, { keyPath: ['scope', 'table', 'key'] });
    opening.result.createObjectStore(GENERATIONS);
  };
  return request(opening).then((db) => {
    const missing = STORES.filter((store) => !db.objectStoreNames.contains(store));
    if (missing.length === 0) return db;
    db.close();
    throw new Error(
      `IndexedDB database "${name}" holds something else (no ${missing.join(', ')} store); use another name`,
    );
  });
}

/**
 * The threads' rows in IndexedDB, the browser's backend (packages/ai/_spec/threads.md § Storage): one database, one
 * record per row keyed [scope, table, key], and one generation per scope. `apply` is one `readwrite` transaction that
 * reads the generation and then only issues puts and deletes from that request's callback, so the transaction never
 * waits on anything but itself (IndexedDB commits a transaction with nothing left to do), and the check and the writes
 * are one step even against a second tab: a tab that wrote the same scope first makes the generation differ and this
 * write is refused as a `RowConflict`, with nothing written. After a commit the scope is told on a `BroadcastChannel`,
 * so other tabs hear of it (`watch`), as the harness's backend does.
 */
export function createIndexedDbRows<T extends string = string>(
  name: string,
): RowStore<T> & { readonly database: string; clear(): Promise<void>; close(): void } {
  const db = open(name);
  const channel = typeof BroadcastChannel === 'undefined' ? undefined : new BroadcastChannel(`threads:${name}`);
  const changed = (scope: string) => channel?.postMessage({ scope });
  return {
    /** Its name: what the kept records' index is fenced by, across instances and tabs (journal.md §1). */
    database: name,
    async load(scope) {
      const tx = (await db).transaction(STORES, 'readonly');
      const [generation, records] = await Promise.all([
        request(tx.objectStore(GENERATIONS).get(scope)) as Promise<number | undefined>,
        request(tx.objectStore(ROWS).getAll(scopeRange(scope))) as Promise<RowRecord[]>,
      ]);
      if (generation === undefined) return undefined;
      return { generation, rows: records.map((r) => ({ table: r.table as T, key: r.key, row: r.row })) };
    },
    async apply(scope, expected, writes: readonly RowWrite<T>[]) {
      const tx = (await db).transaction(STORES, 'readwrite');
      const committed = done(tx);
      let applied = false;
      const generations = tx.objectStore(GENERATIONS);
      const get = generations.get(scope);
      get.onsuccess = () => {
        if (((get.result as number | undefined) ?? 0) !== expected) return;
        applied = true;
        const rows = tx.objectStore(ROWS);
        for (const w of writes) {
          if (w.row) rows.put({ scope, table: w.table, key: w.key, row: w.row } satisfies RowRecord);
          else rows.delete([scope, w.table, w.key]);
        }
        generations.put(expected + 1, scope);
      };
      await committed;
      if (!applied) throw new RowConflict(scope);
      changed(scope);
    },
    async generation(scope) {
      const tx = (await db).transaction(GENERATIONS, 'readonly');
      return ((await request(tx.objectStore(GENERATIONS).get(scope))) as number | undefined) ?? 0;
    },
    async drop(scope) {
      const tx = (await db).transaction(STORES, 'readwrite');
      tx.objectStore(ROWS).delete(scopeRange(scope));
      tx.objectStore(GENERATIONS).delete(scope);
      await done(tx);
      changed(scope);
    },
    watch(listener) {
      if (channel === undefined) return () => undefined;
      const receive = (message: MessageEvent<{ scope: string }>) => listener(message.data.scope);
      channel.addEventListener('message', receive);
      return () => channel.removeEventListener('message', receive);
    },
    /** Forgets every scope. */
    async clear() {
      const tx = (await db).transaction(STORES, 'readwrite');
      tx.objectStore(ROWS).clear();
      tx.objectStore(GENERATIONS).clear();
      await done(tx);
    },
    close() {
      channel?.close();
      void db.then(
        (opened) => opened.close(),
        () => undefined,
      );
    },
  };
}
