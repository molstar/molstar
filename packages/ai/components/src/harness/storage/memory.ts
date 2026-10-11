/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, Message } from '@molstar/ai-agent';
import {
  type BlobRef,
  type Execution,
  HarnessError,
  type HarnessEvent,
  type HarnessId,
  type HarnessObservable,
  type HarnessStorage,
  type HeadRecord,
  type Pointer,
  type RecordStore,
  type ScopeId,
  type Segment,
  type Seq,
  type StepEvent,
  type StorageTx,
} from '@molstar/ai-harness';

/** The Execution is one record per harness, whatever id it is read or written under. */
const EXECUTION_SLOT = 'execution';

export function deepEqual(a: unknown, b: unknown): boolean {
  if (a === b) return true;
  if (typeof a !== 'object' || typeof b !== 'object' || a === null || b === null) return false;
  if (Array.isArray(a) !== Array.isArray(b)) return false;
  const ka = Object.keys(a);
  const kb = Object.keys(b);
  if (ka.length !== kb.length) return false;
  return ka.every(
    (k) => Object.hasOwn(b, k) && deepEqual((a as Record<string, unknown>)[k], (b as Record<string, unknown>)[k]),
  );
}

/** The first index whose event has a seq above `seq` (the list is ascending by seq). */
export function firstAfter(list: readonly HarnessEvent[], seq: number): number {
  let lo = 0;
  let hi = list.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if ((list[mid]?.seq ?? 0) > seq) hi = mid;
    else lo = mid + 1;
  }
  return lo;
}

function push<K>(map: Map<K, HarnessEvent[]>, key: K, event: HarnessEvent): void {
  const list = map.get(key);
  if (list) list.push(event);
  else map.set(key, [event]);
}

/** An append-only event log with its indexes (harness.md §3); each index lists events in seq order. */
export interface Log {
  readonly events: HarnessEvent[];
  readonly steps: Map<string, HarnessEvent[]>;
  readonly turns: Map<string, HarnessEvent[]>;
  readonly keys: Map<string, HarnessEvent[]>;
  readonly compactions: Map<ScopeId, HarnessEvent[]>;
  readonly ended: HarnessEvent[];
}

export const createLog = (): Log => ({
  events: [],
  steps: new Map(),
  turns: new Map(),
  keys: new Map(),
  compactions: new Map(),
  ended: [],
});

export function lastSeq(log: Log): number {
  return log.events.at(-1)?.seq ?? Number.NEGATIVE_INFINITY;
}

export function addEvent(log: Log, event: HarnessEvent): void {
  log.events.push(event);
  if (event.turnId !== undefined) push(log.turns, event.turnId, event);
  if (event.kind === 'step') {
    push(log.steps, event.step, event);
    push(log.keys, event.commitKey, event);
    for (const pointer of new Set((event.compaction ?? []).map((c) => c.pointer))) {
      push(log.compactions, pointer, event);
    }
    if (event.end !== undefined) log.ended.push(event);
  } else {
    if (event.key !== undefined) push(log.keys, event.key, event);
    if (event.summary !== undefined) log.ended.push(event);
  }
}

/** One harness's committed records. */
interface Data {
  readonly execution: Map<string, Execution>;
  readonly messages: Map<string, Message>;
  readonly pointers: Map<string, Pointer>;
  readonly heads: Map<string, HeadRecord>;
  readonly segments: Map<string, Segment>;
  readonly log: Log;
  readonly blobs: Map<BlobRef, Json>;
}

const createData = (): Data => ({
  execution: new Map(),
  messages: new Map(),
  pointers: new Map(),
  heads: new Map(),
  segments: new Map(),
  log: createLog(),
  blobs: new Map(),
});

function recordStore<T>(
  base: Map<string, T>,
  staged: Map<string, T>,
  slot: (id: string) => string,
  immutable: boolean,
): RecordStore<T> {
  return {
    get: async (ids) => ids.map((id) => structuredClone(staged.get(slot(id)) ?? base.get(slot(id)))),
    put(id, value) {
      const key = slot(id);
      const existing = staged.get(key) ?? base.get(key);
      if (immutable && existing !== undefined) {
        if (deepEqual(existing, value)) return;
        throw new HarnessError('integrity', `a stored record cannot change: ${key}`);
      }
      staged.set(key, structuredClone(value));
    },
  };
}

/** The in-memory backend: serializable per harness, buffered writes, structuredClone at the edges. */
export function createMemoryStorage(): HarnessStorage {
  const datas = new Map<HarnessId, Data>();
  const tombstones = new Set<HarnessId>();
  const tails = new Map<HarnessId, Promise<unknown>>();
  const watchers = new Map<HarnessId, Set<{ next: (seq: Seq) => void; error?: (e: unknown) => void }>>();
  const dataOf = (id: HarnessId): Data => {
    const d = datas.get(id) ?? createData();
    datas.set(id, d);
    return d;
  };

  /** Runs `fn` after everything queued before it for the same harness. */
  function serial<T>(id: HarnessId, fn: () => Promise<T>): Promise<T> {
    const result = (tails.get(id) ?? Promise.resolve()).then(fn);
    const tail = result.catch(() => {});
    tails.set(id, tail);
    void tail.then(() => {
      if (tails.get(id) === tail) tails.delete(id);
    });
    return result;
  }

  function emit(id: HarnessId, seq: Seq): void {
    for (const w of [...(watchers.get(id) ?? [])]) {
      try {
        w.next(seq);
      } catch (e) {
        w.error?.(e);
      }
    }
  }

  async function transact<T>(harnessId: HarnessId, run: (tx: StorageTx) => Promise<T>): Promise<T> {
    return serial(harnessId, async () => {
      const base = dataOf(harnessId);
      const staged = createData();
      const stagedTombstones = new Set<HarnessId>();
      const logs = [base.log, staged.log];
      const gather = <E extends HarnessEvent>(
        ordered: Log[],
        limit: number,
        from: (log: Log, n: number) => E[],
      ): E[] => {
        const out: E[] = [];
        for (const log of ordered) if (out.length < limit) out.push(...from(log, limit - out.length));
        return structuredClone(out);
      };
      const events = {
        append(event: HarnessEvent) {
          if (event.seq <= Math.max(lastSeq(base.log), lastSeq(staged.log))) {
            throw new HarnessError('integrity', `event seq ${event.seq} already exists or is out of order`);
          }
          addEvent(staged.log, structuredClone(event));
        },
        after: async (seq: Seq, limit: number) =>
          gather(logs, limit, (l, n) => l.events.slice(firstAfter(l.events, seq)).slice(0, n)),
        byStep: async (ids: readonly string[]) =>
          gather(
            logs,
            Number.POSITIVE_INFINITY,
            (l) => [...new Set(ids)].flatMap((id) => l.steps.get(id) ?? []) as StepEvent[],
          ).sort((a, b) => a.seq - b.seq),
        byTurn: async (turnId: string, after: Seq | undefined, limit: number) =>
          gather(logs, limit, (l, n) => {
            const list = l.turns.get(turnId) ?? [];
            return list.slice(firstAfter(list, after ?? Number.NEGATIVE_INFINITY)).slice(0, n);
          }),
        endedTurns: async (before: Seq | undefined, limit: number) =>
          gather([...logs].reverse(), limit, (l, n) => {
            const end = before === undefined ? l.ended.length : firstAfter(l.ended, before - 1);
            return l.ended.slice(Math.max(0, end - n), end).reverse();
          }),
        byKey: async (key: string) => gather(logs, Number.POSITIVE_INFINITY, (l) => l.keys.get(key) ?? []),
        lastCompaction: async (pointer: ScopeId) => {
          const last = staged.log.compactions.get(pointer)?.at(-1) ?? base.log.compactions.get(pointer)?.at(-1);
          return structuredClone(last as StepEvent | undefined);
        },
      };
      const same = (id: string) => id;
      const tx: StorageTx = {
        execution: recordStore(base.execution, staged.execution, () => EXECUTION_SLOT, false),
        events,
        messages: recordStore(base.messages, staged.messages, same, true),
        pointers: recordStore(base.pointers, staged.pointers, same, false),
        heads: recordStore(base.heads, staged.heads, same, true),
        segments: recordStore(base.segments, staged.segments, same, true),
        store: {
          tombstone: async (id) => tombstones.has(id) || stagedTombstones.has(id),
          setTombstone: (id) => void stagedTombstones.add(id),
        },
      };
      const out = await run(tx);
      const before = base.execution.get(EXECUTION_SLOT)?.seq ?? 0;
      for (const key of ['execution', 'messages', 'pointers', 'heads', 'segments'] as const) {
        for (const [k, v] of staged[key]) (base[key] as Map<string, unknown>).set(k, v);
      }
      for (const event of staged.log.events) addEvent(base.log, event);
      for (const id of stagedTombstones) tombstones.add(id);
      const after = base.execution.get(EXECUTION_SLOT)?.seq ?? 0;
      if (after !== before) emit(harnessId, after);
      return out;
    });
  }

  return {
    transact,
    async putBlobs(harnessId, blobs) {
      for (const { ref, body } of blobs) dataOf(harnessId).blobs.set(ref, structuredClone(body));
    },
    async getBlobs(harnessId, refs) {
      return refs.map((ref) => structuredClone(dataOf(harnessId).blobs.get(ref)));
    },
    watch(harnessId): HarnessObservable<Seq> {
      return {
        subscribe(next, error) {
          const w = error ? { next, error } : { next };
          const set = watchers.get(harnessId) ?? new Set();
          watchers.set(harnessId, set);
          set.add(w);
          try {
            next(dataOf(harnessId).execution.get(EXECUTION_SLOT)?.seq ?? 0);
          } catch (e) {
            error?.(e);
          }
          return () => void set.delete(w);
        },
      };
    },
    drop: (harnessId) => serial(harnessId, async () => void datas.delete(harnessId)),
  };
}
