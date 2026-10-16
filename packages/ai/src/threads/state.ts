/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ThreadDeltaBatch, ThreadEntry } from './domain.js';
import { threadCommandId, threadEntryId, threadId, threadOwnerId } from './ids.js';
import type { ThreadKernelThread, ThreadKernelReceipt, ThreadKernelStream } from './kernel.js';
import type { ThreadTable, Row, RowWrite } from './storage.js';

/** A stream as the host keeps it: the kernel's view, where its deltas are, when a preview expires. */
export type Stream = ThreadKernelStream & { readonly streamId: string; readonly expiresAt?: number };

/**
 * One thread as loaded: what anyfold's in-memory host keeps for all threads in one state, here per
 * thread, so a thread is its own scope of rows, loaded when first used. Maps are replaced, never
 * changed, by a command (copy on write): a write set is the values that differ by reference.
 */
export interface ThreadState {
  /** The scope's generation as last read or written. */
  generation: number;
  thread: ThreadKernelThread | undefined;
  /** By entry id. */
  entries: Map<string, ThreadEntry>;
  /** Entry ids by position. */
  order: string[];
  /** By entry id. */
  streams: Map<string, Stream>;
  /** By stream id. */
  deltas: Map<string, ThreadDeltaBatch[]>;
  /** By command id. */
  receipts: Map<string, ThreadKernelReceipt>;
}

export function emptyState(): ThreadState {
  return {
    generation: 0,
    thread: undefined,
    entries: new Map(),
    order: [],
    streams: new Map(),
    deltas: new Map(),
    receipts: new Map(),
  };
}

/** A draft for one command: new containers over the same values. */
export function draftOf(state: ThreadState): ThreadState {
  return {
    generation: state.generation,
    thread: state.thread,
    entries: new Map(state.entries),
    order: [...state.order],
    streams: new Map(state.streams),
    deltas: new Map(state.deltas),
    receipts: new Map(state.receipts),
  };
}

const THREAD_KEY = 'thread';

/** A thread's state from its scope's rows (JSON restores the branded ids as plain strings). */
export function stateOf(
  loaded: { generation: number; rows: readonly RowWrite<ThreadTable>[] } | undefined,
): ThreadState {
  const state = emptyState();
  if (!loaded) return state;
  state.generation = loaded.generation;
  const positions: [number, string][] = [];
  for (const { table, key, row } of loaded.rows) {
    if (!row) continue;
    const value = JSON.parse(String(row.json));
    switch (table) {
      case 'thread':
        state.thread = { ...value, threadId: threadId(value.threadId), ownerId: threadOwnerId(value.ownerId) };
        break;
      case 'entry':
        state.entries.set(key, { ...value, entryId: threadEntryId(value.entryId), threadId: threadId(value.threadId) });
        positions.push([value.position, key]);
        break;
      case 'stream':
        state.streams.set(key, { ...value, entryId: threadEntryId(value.entryId) });
        break;
      case 'deltas':
        state.deltas.set(
          key,
          (value as ThreadDeltaBatch[]).map((b) => ({ ...b, entryId: threadEntryId(b.entryId) })),
        );
        break;
      case 'receipt':
        state.receipts.set(key, { ...value, commandId: threadCommandId(value.commandId) });
        break;
    }
  }
  state.order = positions.sort((a, b) => a[0] - b[0]).map(([, id]) => id);
  return state;
}

const row = (value: unknown): Row => ({ json: JSON.stringify(value) });

/** What changed from `before` to `after`, as a write set. */
export function writesOf(before: ThreadState, after: ThreadState): RowWrite<ThreadTable>[] {
  const writes: RowWrite<ThreadTable>[] = [];
  if (after.thread !== before.thread)
    writes.push({ table: 'thread', key: THREAD_KEY, row: after.thread ? row(after.thread) : undefined });
  const table = <V>(name: ThreadTable, a: Map<string, V>, b: Map<string, V>) => {
    for (const [key, value] of b) if (a.get(key) !== value) writes.push({ table: name, key, row: row(value) });
    for (const key of a.keys()) if (!b.has(key)) writes.push({ table: name, key, row: undefined });
  };
  table('entry', before.entries, after.entries);
  table('stream', before.streams, after.streams);
  table('deltas', before.deltas, after.deltas);
  table('receipt', before.receipts, after.receipts);
  return writes;
}
