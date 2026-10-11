/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { FiberDeltaBatch, FiberEntry } from './domain.js';
import { fiberCommandId, fiberEntryId, fiberId, fiberOwnerId } from './ids.js';
import type { FiberKernelFiber, FiberKernelReceipt, FiberKernelStream } from './kernel.js';
import type { FiberTable, Row, RowWrite } from './storage.js';

/** A stream as the host keeps it: the kernel's view, where its deltas are, when a preview expires. */
export type Stream = FiberKernelStream & { readonly streamId: string; readonly expiresAt?: number };

/**
 * One fiber as loaded: what anyfold's in-memory host keeps for all fibers in one state, here per
 * fiber, so a fiber is its own scope of rows, loaded when first used. Maps are replaced, never
 * changed, by a command (copy on write): a write set is the values that differ by reference.
 */
export interface FiberState {
  /** The scope's generation as last read or written. */
  generation: number;
  fiber: FiberKernelFiber | undefined;
  /** By entry id. */
  entries: Map<string, FiberEntry>;
  /** Entry ids by position. */
  order: string[];
  /** By entry id. */
  streams: Map<string, Stream>;
  /** By stream id. */
  deltas: Map<string, FiberDeltaBatch[]>;
  /** By command id. */
  receipts: Map<string, FiberKernelReceipt>;
}

export function emptyState(): FiberState {
  return {
    generation: 0,
    fiber: undefined,
    entries: new Map(),
    order: [],
    streams: new Map(),
    deltas: new Map(),
    receipts: new Map(),
  };
}

/** A draft for one command: new containers over the same values. */
export function draftOf(state: FiberState): FiberState {
  return {
    generation: state.generation,
    fiber: state.fiber,
    entries: new Map(state.entries),
    order: [...state.order],
    streams: new Map(state.streams),
    deltas: new Map(state.deltas),
    receipts: new Map(state.receipts),
  };
}

const FIBER_KEY = 'fiber';

/** A fiber's state from its scope's rows (JSON restores the branded ids as plain strings). */
export function stateOf(loaded: { generation: number; rows: readonly RowWrite<FiberTable>[] } | undefined): FiberState {
  const state = emptyState();
  if (!loaded) return state;
  state.generation = loaded.generation;
  const positions: [number, string][] = [];
  for (const { table, key, row } of loaded.rows) {
    if (!row) continue;
    const value = JSON.parse(String(row.json));
    switch (table) {
      case 'fiber':
        state.fiber = { ...value, fiberId: fiberId(value.fiberId), ownerId: fiberOwnerId(value.ownerId) };
        break;
      case 'entry':
        state.entries.set(key, { ...value, entryId: fiberEntryId(value.entryId), fiberId: fiberId(value.fiberId) });
        positions.push([value.position, key]);
        break;
      case 'stream':
        state.streams.set(key, { ...value, entryId: fiberEntryId(value.entryId) });
        break;
      case 'deltas':
        state.deltas.set(
          key,
          (value as FiberDeltaBatch[]).map((b) => ({ ...b, entryId: fiberEntryId(b.entryId) })),
        );
        break;
      case 'receipt':
        state.receipts.set(key, { ...value, commandId: fiberCommandId(value.commandId) });
        break;
    }
  }
  state.order = positions.sort((a, b) => a[0] - b[0]).map(([, id]) => id);
  return state;
}

const row = (value: unknown): Row => ({ json: JSON.stringify(value) });

/** What changed from `before` to `after`, as a write set. */
export function writesOf(before: FiberState, after: FiberState): RowWrite<FiberTable>[] {
  const writes: RowWrite<FiberTable>[] = [];
  if (after.fiber !== before.fiber)
    writes.push({ table: 'fiber', key: FIBER_KEY, row: after.fiber ? row(after.fiber) : undefined });
  const table = <V>(name: FiberTable, a: Map<string, V>, b: Map<string, V>) => {
    for (const [key, value] of b) if (a.get(key) !== value) writes.push({ table: name, key, row: row(value) });
    for (const key of a.keys()) if (!b.has(key)) writes.push({ table: name, key, row: undefined });
  };
  table('entry', before.entries, after.entries);
  table('stream', before.streams, after.streams);
  table('deltas', before.deltas, after.deltas);
  table('receipt', before.receipts, after.receipts);
  return writes;
}
