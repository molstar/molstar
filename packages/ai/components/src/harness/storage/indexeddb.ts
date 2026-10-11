/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';
import type { HarnessEvent, HarnessStorage, StepEvent } from '@molstar/ai-harness';
import { createDurableStorage, type DurableBackend } from './durable.js';

/** An event as stored, with the fields its indexes use (absent fields are not indexed). */
interface StoredEvent {
  readonly h: string;
  readonly seq: number;
  readonly event: HarnessEvent;
  readonly step?: string;
  readonly turn?: string;
  readonly key?: string;
  readonly ended?: 1;
  /** `[h, pointer, seq]` per compacted pointer (a multi-entry index). */
  readonly compactions?: readonly [string, string, number][];
}

const STORES = ['meta', 'records', 'events', 'blobs'] as const;
const MAX = Number.MAX_SAFE_INTEGER;
/** `getAll` counts are unsigned 32-bit; anything larger means all. */
const count = (limit: number) => (limit >= 0xffffffff ? undefined : limit);

const request = <T>(r: IDBRequest<T>): Promise<T> =>
  new Promise((resolve, reject) => {
    r.onsuccess = () => resolve(r.result);
    r.onerror = () => reject(r.error);
  });

function open(name: string): Promise<IDBDatabase> {
  const opening = indexedDB.open(name, 1);
  opening.onupgradeneeded = () => {
    const db = opening.result;
    db.createObjectStore('meta');
    db.createObjectStore('records');
    db.createObjectStore('blobs');
    const events = db.createObjectStore('events', { keyPath: ['h', 'seq'] });
    events.createIndex('step', ['h', 'step']);
    events.createIndex('turn', ['h', 'turn', 'seq']);
    events.createIndex('key', ['h', 'key']);
    events.createIndex('ended', ['h', 'ended', 'seq']);
    events.createIndex('compactions', 'compactions', { multiEntry: true });
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

function stored(h: string, event: HarnessEvent): StoredEvent {
  const step = event.kind === 'step' ? event : undefined;
  const key = event.kind === 'step' ? event.commitKey : event.key;
  const ended = event.kind === 'step' ? event.end !== undefined : event.summary !== undefined;
  const pointers = [...new Set((step?.compaction ?? []).map((c) => c.pointer))];
  return {
    h,
    seq: event.seq,
    event,
    ...(step !== undefined ? { step: step.step } : {}),
    ...(event.turnId !== undefined ? { turn: event.turnId } : {}),
    ...(key !== undefined ? { key } : {}),
    ...(ended ? { ended: 1 as const } : {}),
    ...(pointers.length > 0 ? { compactions: pointers.map((p): [string, string, number] => [h, p, event.seq]) } : {}),
  };
}

/**
 * The browser's backend (harness.md §14): one IndexedDB database. Writes apply in one `readwrite`
 * transaction that first checks the harness's version, issued without awaiting anything else, so
 * the transaction never closes early; other tabs hear of them on a `BroadcastChannel`.
 */
export function createIndexedDbBackend(name: string): DurableBackend & { close(): void } {
  const db = open(name);
  const channel = typeof BroadcastChannel === 'undefined' ? undefined : new BroadcastChannel(`harness:${name}`);
  const read = async <T>(store: (typeof STORES)[number], use: (s: IDBObjectStore) => IDBRequest<T>) =>
    await request(use((await db).transaction(store, 'readonly').objectStore(store)));
  const range = (
    index: string,
    lower: unknown[],
    upper: unknown[],
    limit: number,
    direction: IDBCursorDirection = 'next',
  ) =>
    read('events', (s) => s.index(index).getAll(IDBKeyRange.bound(lower, upper), count(limit))).then((rows) => {
      const events = (rows as StoredEvent[]).map((row) => row.event);
      return direction === 'prev' ? events.reverse() : events;
    });

  return {
    async get(h, kind, ids) {
      const store = (await db).transaction('records', 'readonly').objectStore('records');
      return await Promise.all(ids.map((id) => request(store.get([h, kind, id]))));
    },
    after: async (h, seq, limit) =>
      (
        (await read('events', (s) =>
          s.getAll(IDBKeyRange.bound([h, seq], [h, MAX], true), count(limit)),
        )) as StoredEvent[]
      ).map((row) => row.event),
    async byStep(h, ids) {
      const found = await Promise.all([...new Set(ids)].map((id) => range('step', [h, id], [h, id], MAX)));
      return (found.flat() as StepEvent[]).sort((a, b) => a.seq - b.seq);
    },
    byTurn: (h, turn, after, limit) => range('turn', [h, turn, (after ?? -1) + 1], [h, turn, MAX], limit),
    async endedTurns(h, before, limit) {
      const rows: HarnessEvent[] = [];
      const index = (await db).transaction('events', 'readonly').objectStore('events').index('ended');
      const cursor = index.openCursor(IDBKeyRange.bound([h, 1, -1], [h, 1, (before ?? MAX) - 1]), 'prev');
      return await new Promise((resolve, reject) => {
        cursor.onerror = () => reject(cursor.error);
        cursor.onsuccess = () => {
          const at = cursor.result;
          if (at === null || rows.length >= limit) return resolve(rows);
          rows.push((at.value as StoredEvent).event);
          at.continue();
        };
      });
    },
    byKey: (h, key) => range('key', [h, key], [h, key], MAX),
    async lastCompaction(h, pointer) {
      const index = (await db).transaction('events', 'readonly').objectStore('events').index('compactions');
      const cursor = await request(index.openCursor(IDBKeyRange.bound([h, pointer, -1], [h, pointer, MAX]), 'prev'));
      return (cursor?.value as StoredEvent | undefined)?.event as StepEvent | undefined;
    },
    tombstone: async (h) =>
      ((await read('meta', (s) => s.get(h))) as { tombstone?: boolean } | undefined)?.tombstone === true,
    version: async (h) => ((await read('meta', (s) => s.get(h))) as { version?: number } | undefined)?.version ?? 0,
    async apply(h, expected, writes) {
      const tx = (await db).transaction(STORES, 'readwrite');
      const done = new Promise<void>((resolve, reject) => {
        tx.oncomplete = () => resolve();
        tx.onerror = () => reject(tx.error);
        tx.onabort = () => reject(tx.error);
      });
      let applied = false;
      const meta = tx.objectStore('meta').get(h);
      meta.onsuccess = () => {
        const current = (meta.result as { version: number; tombstone?: boolean } | undefined) ?? { version: 0 };
        if (current.version !== expected) return;
        applied = true;
        const tombstone = current.tombstone === true || writes.tombstone;
        tx.objectStore('meta').put({ version: expected + 1, ...(tombstone ? { tombstone } : {}) }, h);
        for (const { kind, id, value } of writes.records) tx.objectStore('records').put(value, [h, kind, id]);
        for (const event of writes.events) tx.objectStore('events').add(stored(h, event));
      };
      await done;
      const execution = writes.records.find((record) => record.kind === 'execution')?.value as
        | { seq: number }
        | undefined;
      if (applied && execution !== undefined) channel?.postMessage({ h, seq: execution.seq });
      return applied;
    },
    async putBlobs(h, blobs) {
      const tx = (await db).transaction('blobs', 'readwrite');
      for (const { ref, body } of blobs) tx.objectStore('blobs').put(body, [h, ref]);
      await new Promise<void>((resolve, reject) => {
        tx.oncomplete = () => resolve();
        tx.onerror = () => reject(tx.error);
      });
    },
    async getBlobs(h, refs) {
      const store = (await db).transaction('blobs', 'readonly').objectStore('blobs');
      return await Promise.all(refs.map((ref) => request(store.get([h, ref])) as Promise<Json | undefined>));
    },
    async drop(h) {
      const tx = (await db).transaction(['records', 'events', 'blobs'], 'readwrite');
      tx.objectStore('records').delete(IDBKeyRange.bound([h], [h, []]));
      tx.objectStore('events').delete(IDBKeyRange.bound([h, -1], [h, MAX]));
      tx.objectStore('blobs').delete(IDBKeyRange.bound([h], [h, []]));
      await new Promise<void>((resolve, reject) => {
        tx.oncomplete = () => resolve();
        tx.onerror = () => reject(tx.error);
      });
    },
    onChange(listener) {
      if (channel === undefined) return () => undefined;
      const receive = (message: MessageEvent<{ h: string; seq: number }>) => listener(message.data.h, message.data.seq);
      channel.addEventListener('message', receive);
      return () => channel.removeEventListener('message', receive);
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

/** A durable storage in one IndexedDB database. */
export function createIndexedDbStorage(name: string): HarnessStorage & { close(): void } {
  const backend = createIndexedDbBackend(name);
  return Object.assign(createDurableStorage(backend), { close: () => backend.close() });
}
