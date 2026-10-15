/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { canonicalJson, type Json } from '@molstar/ai-agent';
import type { EventLog, HarnessStorage, RecordStore, StorageTx } from '@molstar/ai-harness';

export type RecordKind = 'execution' | 'events' | 'messages' | 'pointers' | 'heads' | 'segments' | 'blobs' | 'store';

export interface StorageCounts {
  transactions: number;
  reads: Record<RecordKind, number>;
  writes: Record<RecordKind, number>;
  bytesRead: number;
  bytesWritten: number;
}

const KINDS: readonly RecordKind[] = [
  'execution',
  'events',
  'messages',
  'pointers',
  'heads',
  'segments',
  'blobs',
  'store',
];

export function emptyCounts(): StorageCounts {
  const zero = () => Object.fromEntries(KINDS.map((kind) => [kind, 0])) as Record<RecordKind, number>;
  return { transactions: 0, reads: zero(), writes: zero(), bytesRead: 0, bytesWritten: 0 };
}

/** `after` minus `before`, field by field. */
export function difference(after: StorageCounts, before: StorageCounts): StorageCounts {
  const by = (a: Record<RecordKind, number>, b: Record<RecordKind, number>) =>
    Object.fromEntries(KINDS.map((kind) => [kind, a[kind] - b[kind]])) as Record<RecordKind, number>;
  return {
    transactions: after.transactions - before.transactions,
    reads: by(after.reads, before.reads),
    writes: by(after.writes, before.writes),
    bytesRead: after.bytesRead - before.bytesRead,
    bytesWritten: after.bytesWritten - before.bytesWritten,
  };
}

export function totalReads(counts: StorageCounts): number {
  return KINDS.reduce((sum, kind) => sum + counts.reads[kind], 0);
}
export function totalWrites(counts: StorageCounts): number {
  return KINDS.reduce((sum, kind) => sum + counts.writes[kind], 0);
}

const size = (value: unknown): number => (value === undefined ? 0 : canonicalJson(value as Json).length);

/** Wraps a backend so every transaction, record read and write and their bytes are counted. */
export function meterStorage(storage: HarnessStorage): { storage: HarnessStorage; counts(): StorageCounts } {
  const counts = emptyCounts();
  const read = (kind: RecordKind, values: readonly unknown[]) => {
    for (const value of values) {
      if (value === undefined) continue;
      counts.reads[kind] += 1;
      counts.bytesRead += size(value);
    }
  };
  const write = (kind: RecordKind, value: unknown) => {
    counts.writes[kind] += 1;
    counts.bytesWritten += size(value);
  };
  const store = <T>(kind: RecordKind, inner: RecordStore<T>): RecordStore<T> => ({
    async get(ids) {
      const values = await inner.get(ids);
      read(kind, values);
      return values;
    },
    put(id, value) {
      write(kind, value);
      inner.put(id, value);
    },
  });
  const events = (inner: EventLog): EventLog => ({
    append(event) {
      write('events', event);
      inner.append(event);
    },
    after: async (seq, limit) => {
      const out = await inner.after(seq, limit);
      read('events', out);
      return out;
    },
    byStep: async (ids) => {
      const out = await inner.byStep(ids);
      read('events', out);
      return out;
    },
    byTurn: async (turnId, after, limit) => {
      const out = await inner.byTurn(turnId, after, limit);
      read('events', out);
      return out;
    },
    endedTurns: async (before, limit) => {
      const out = await inner.endedTurns(before, limit);
      read('events', out);
      return out;
    },
    byKey: async (key) => {
      const out = await inner.byKey(key);
      read('events', out);
      return out;
    },
    lastCompaction: async (pointer) => {
      const out = await inner.lastCompaction(pointer);
      read('events', out === undefined ? [] : [out]);
      return out;
    },
  });
  const metered: HarnessStorage = {
    transact(harnessId, run) {
      counts.transactions += 1;
      return storage.transact(harnessId, (tx) =>
        run({
          execution: store('execution', tx.execution),
          events: events(tx.events),
          messages: store('messages', tx.messages),
          pointers: store('pointers', tx.pointers),
          heads: store('heads', tx.heads),
          segments: store('segments', tx.segments),
          store: {
            tombstone: async (id) => {
              counts.reads.store += 1;
              return await tx.store.tombstone(id);
            },
            setTombstone: (id) => {
              counts.writes.store += 1;
              tx.store.setTombstone(id);
            },
          },
        } satisfies StorageTx),
      );
    },
    async putBlobs(harnessId, blobs) {
      for (const blob of blobs) write('blobs', blob.body);
      await storage.putBlobs(harnessId, blobs);
    },
    async getBlobs(harnessId, refs) {
      const out = await storage.getBlobs(harnessId, refs);
      read('blobs', out);
      return out;
    },
    watch: (harnessId) => storage.watch(harnessId),
    drop: (harnessId) => storage.drop(harnessId),
  };
  return { storage: metered, counts: () => structuredClone(counts) };
}
