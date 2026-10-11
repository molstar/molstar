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
  type RecordStore,
  type ScopeId,
  type Seq,
  type StepEvent,
  type StepId,
  type StorageTx,
  type TurnId,
} from '@molstar/ai-harness';
import { addEvent, createLog, deepEqual, firstAfter, type Log, lastSeq } from './memory.js';

export type RecordKind = 'execution' | 'messages' | 'pointers' | 'heads' | 'segments';
const IMMUTABLE: ReadonlySet<RecordKind> = new Set(['messages', 'heads', 'segments']);

/** One transaction's writes, applied all at once. The Execution's id is the harness id. */
export interface WriteSet {
  readonly records: readonly { readonly kind: RecordKind; readonly id: string; readonly value: unknown }[];
  readonly events: readonly HarnessEvent[];
  readonly tombstone: boolean;
}

/**
 * What a durable backend provides: committed reads, the event log's indexed queries, and one atomic
 * apply guarded by the harness's write version. The transaction logic is `createDurableStorage`'s.
 */
export interface DurableBackend {
  get(harnessId: HarnessId, kind: RecordKind, ids: readonly string[]): Promise<(unknown | undefined)[]>;
  after(harnessId: HarnessId, seq: Seq, limit: number): Promise<HarnessEvent[]>;
  byStep(harnessId: HarnessId, ids: readonly StepId[]): Promise<StepEvent[]>;
  byTurn(harnessId: HarnessId, turnId: TurnId, after: Seq | undefined, limit: number): Promise<HarnessEvent[]>;
  /** Newest first, below `before`. */
  endedTurns(harnessId: HarnessId, before: Seq | undefined, limit: number): Promise<HarnessEvent[]>;
  byKey(harnessId: HarnessId, key: string): Promise<HarnessEvent[]>;
  lastCompaction(harnessId: HarnessId, pointer: ScopeId): Promise<StepEvent | undefined>;
  tombstone(harnessId: HarnessId): Promise<boolean>;
  /** The harness's write version: 0 before its first write. */
  version(harnessId: HarnessId): Promise<number>;
  /** Applies the writes and bumps the version if it is still `expected`; false otherwise. */
  apply(harnessId: HarnessId, expected: number, writes: WriteSet): Promise<boolean>;
  putBlobs(harnessId: HarnessId, blobs: readonly { readonly ref: BlobRef; readonly body: Json }[]): Promise<void>;
  getBlobs(harnessId: HarnessId, refs: readonly BlobRef[]): Promise<(Json | undefined)[]>;
  drop(harnessId: HarnessId): Promise<void>;
  /** Changes applied by other processes (another tab), if the backend can see them. */
  onChange?(listener: (harnessId: HarnessId, seq: Seq) => void): () => void;
}

/** Events that indexes would return, from the staged log, after (or, for `endedTurns`, before) the backend's. */
function stagedEvents(log: Log) {
  return {
    after: (seq: Seq) => log.events.slice(firstAfter(log.events, seq)),
    byStep: (ids: readonly string[]) => [...new Set(ids)].flatMap((id) => log.steps.get(id) ?? []) as StepEvent[],
    byTurn: (turnId: string, after: Seq | undefined) => {
      const list = log.turns.get(turnId) ?? [];
      return list.slice(firstAfter(list, after ?? Number.NEGATIVE_INFINITY));
    },
    ended: (before: Seq | undefined) =>
      log.ended.filter((event) => before === undefined || event.seq < before).reverse(),
  };
}

/**
 * A durable `HarnessStorage` over a backend (harness.md §14): each transaction reads through its own staged
 * writes, then applies them in one atomic write guarded by the harness's version, and runs again on a
 * conflict (another tab or process wrote in between). Within this process, transactions are serial
 * per harness, so conflicts come only from outside it.
 */
export function createDurableStorage(backend: DurableBackend): HarnessStorage {
  const tails = new Map<HarnessId, Promise<unknown>>();
  const watchers = new Map<HarnessId, Set<{ next(seq: Seq): void; error?(e: unknown): void }>>();

  const emit = (harnessId: HarnessId, seq: Seq) => {
    for (const watcher of [...(watchers.get(harnessId) ?? [])]) {
      try {
        watcher.next(seq);
      } catch (error) {
        watcher.error?.(error);
      }
    }
  };
  backend.onChange?.(emit);

  function serial<T>(harnessId: HarnessId, work: () => Promise<T>): Promise<T> {
    const result = (tails.get(harnessId) ?? Promise.resolve()).then(work);
    const tail = result.catch(() => undefined);
    tails.set(harnessId, tail);
    void tail.then(() => {
      if (tails.get(harnessId) === tail) tails.delete(harnessId);
    });
    return result;
  }

  async function attempt<T>(harnessId: HarnessId, run: (tx: StorageTx) => Promise<T>) {
    const version = await backend.version(harnessId);
    const staged = new Map<string, { kind: RecordKind; id: string; value: unknown }>();
    const log = createLog();
    let tombstone = false;
    const store = <V>(kind: RecordKind): RecordStore<V> => ({
      async get(ids) {
        const keys = ids.map((id) => `${kind}:${kind === 'execution' ? harnessId : id}`);
        const missing = ids.filter((_, i) => !staged.has(keys[i] as string));
        const found = await backend.get(harnessId, kind, kind === 'execution' ? missing.map(() => harnessId) : missing);
        let next = 0;
        return keys.map((key) =>
          structuredClone((staged.has(key) ? staged.get(key)?.value : found[next++]) as V | undefined),
        );
      },
      put(id, value) {
        const key = `${kind}:${kind === 'execution' ? harnessId : id}`;
        const prior = staged.get(key);
        if (IMMUTABLE.has(kind) && prior !== undefined && !deepEqual(prior.value, value))
          throw new HarnessError('integrity', `a stored record cannot change: ${key}`);
        staged.set(key, { kind, id: kind === 'execution' ? harnessId : id, value: structuredClone(value) });
      },
    });
    const local = stagedEvents(log);
    const tx: StorageTx = {
      execution: store<Execution>('execution'),
      messages: store<Message>('messages'),
      pointers: store('pointers'),
      heads: store('heads'),
      segments: store('segments'),
      events: {
        append(event) {
          if (log.events.length > 0 && event.seq <= lastSeq(log))
            throw new HarnessError('integrity', `event seq ${event.seq} is out of order`);
          addEvent(log, structuredClone(event));
        },
        after: async (seq, limit) =>
          [...(await backend.after(harnessId, seq, limit)), ...local.after(seq)].slice(0, limit),
        byStep: async (ids) => [...(await backend.byStep(harnessId, ids)), ...local.byStep(ids)],
        byTurn: async (turnId, after, limit) =>
          [...(await backend.byTurn(harnessId, turnId, after, limit)), ...local.byTurn(turnId, after)].slice(0, limit),
        endedTurns: async (before, limit) => {
          const mine = local.ended(before);
          return [...mine, ...(await backend.endedTurns(harnessId, before, Math.max(0, limit - mine.length)))].slice(
            0,
            limit,
          );
        },
        byKey: async (key) => [...(await backend.byKey(harnessId, key)), ...(log.keys.get(key) ?? [])],
        lastCompaction: async (pointer) =>
          (log.compactions.get(pointer)?.at(-1) as StepEvent | undefined) ??
          (await backend.lastCompaction(harnessId, pointer)),
      },
      store: {
        tombstone: async (id) => (id === harnessId && tombstone) || (await backend.tombstone(id)),
        setTombstone: (id) => {
          if (id !== harnessId) throw new HarnessError('invalid', 'a transaction writes only its own harness');
          tombstone = true;
        },
      },
    };
    const out = await run(tx);
    if (staged.size === 0 && log.events.length === 0 && !tombstone) return { out };
    const records = [...staged.values()];
    for (const kind of IMMUTABLE) {
      const writes = records.filter((record) => record.kind === kind);
      const existing = await backend.get(
        harnessId,
        kind,
        writes.map((record) => record.id),
      );
      writes.forEach((record, i) => {
        if (existing[i] !== undefined && !deepEqual(existing[i], record.value))
          throw new HarnessError('integrity', `a stored record cannot change: ${kind}:${record.id}`);
      });
    }
    const applied = await backend.apply(harnessId, version, { records, events: log.events, tombstone });
    if (!applied) return undefined;
    const execution = staged.get(`execution:${harnessId}`)?.value as Execution | undefined;
    if (execution !== undefined) emit(harnessId, execution.seq);
    return { out };
  }

  return {
    transact: (harnessId, run) =>
      serial(harnessId, async () => {
        for (let tries = 1; ; tries += 1) {
          const done = await attempt(harnessId, run);
          if (done !== undefined) return done.out;
          if (tries >= 10) throw new HarnessError('unknown-outcome', 'the transaction kept conflicting');
        }
      }),
    putBlobs: (harnessId, blobs) => backend.putBlobs(harnessId, blobs),
    getBlobs: (harnessId, refs) => backend.getBlobs(harnessId, refs),
    drop: (harnessId) => backend.drop(harnessId),
    watch(harnessId): HarnessObservable<Seq> {
      return {
        subscribe(next, error) {
          const watcher = error !== undefined ? { next, error } : { next };
          const set = watchers.get(harnessId) ?? new Set();
          watchers.set(harnessId, set);
          set.add(watcher);
          backend.get(harnessId, 'execution', [harnessId]).then(
            ([execution]) => {
              if (set.has(watcher)) next((execution as Execution | undefined)?.seq ?? 0);
            },
            (e: unknown) => error?.(e),
          );
          return () => void set.delete(watcher);
        },
      };
    },
  };
}
