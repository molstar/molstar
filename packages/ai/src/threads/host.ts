/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Clock } from '../messages/index.js';
import { cloneData } from './data.js';
import type {
  ThreadCommandResult,
  ThreadDeltaBatch,
  ThreadEntry,
  ThreadFailure,
  ThreadPort,
  ThreadRemoveResult,
} from './domain.js';
import { type ThreadId, threadId } from './ids.js';
import {
  type ThreadKernelThread,
  type ThreadKernelReceipt,
  type ThreadKernelStream,
  threadCommandFingerprint,
  threadPageTakeBound,
  threadStoredEnvelopeBytes,
  threadWindowRangeStart,
  prepareThreadAbort,
  prepareThreadAppendDelta,
  prepareThreadAppendSettled,
  prepareThreadBegin,
  prepareThreadCreate,
  prepareThreadHeartbeat,
  prepareThreadInactivityAbort,
  prepareThreadRemove,
  prepareThreadSettle,
  projectThreadDeltas,
  projectThreadWindow,
} from './kernel.js';
import { createThreadLimits } from './limits.js';
import { draftOf, type ThreadState, type Stream, stateOf, writesOf } from './state.js';
import { type ThreadRows, RowConflict, type Unsubscribe } from './storage.js';

/**
 * The threads' host (threads.md), written from anyfold's in-memory host and its Convex host as
 * Stories ported it: the kernel decides every command; the host keeps each thread's state and writes what a command
 * changed through to its rows before the command returns (so an accepted command is durable). Unlike the in-memory
 * host, a thread is its own scope of rows, loaded when first used and kept in memory after.
 *
 * Commands run one at a time. A replayed command id is found in its thread's receipts, and in the receipts of the
 * other threads loaded here (a command id reused for another thread is a conflict); callers name commands after
 * their thread, so ids never meet across threads.
 */
export interface ThreadHost extends ThreadPort {
  /** Hears which thread a command of this host changed. */
  subscribe(listener: (threadId: ThreadId) => void): Unsubscribe;
  /**
   * Hears which thread was changed through another handle to the same rows (another tab, another process), when
   * the rows can tell (`RowStore.watch`); otherwise never. The next command or read on that thread sees the change
   * with or without this.
   */
  watch(listener: (threadId: ThreadId) => void): Unsubscribe;
  /** Forgets a thread entirely: its rows and what is loaded (its owner was removed). */
  forget(threadId: ThreadId): Promise<void>;
}

export interface ThreadHostOptions {
  rows: ThreadRows | Promise<ThreadRows>;
  /** Only `now` is read. */
  clock?: Pick<Clock, 'now'>;
  /** A thread's scope in the rows; default `thread:<id>`. */
  scopeOf?: (threadId: ThreadId) => string;
  /**
   * The thread of a scope the rows report as changed (`watch`); default: the threads this host has used, and with
   * the default `scopeOf` any `thread:<id>`.
   */
  threadOf?: (scope: string) => ThreadId | undefined;
}

const MAX_CONFLICTS = 8;
/** A read of a thread this store does not hold. */
const MISSING: ThreadFailure = { status: 'not-found', reason: 'thread-missing' };
const DEFAULT_SCOPE_PREFIX = 'thread:';

export function createThreadHost(options: ThreadHostOptions): ThreadHost {
  const clock = options.clock ?? { now: () => Date.now() };
  const defaultScope = options.scopeOf === undefined;
  const known = new Map<string, ThreadId>();
  const scopeOf = (id: ThreadId): string => {
    const scope = options.scopeOf?.(id) ?? `${DEFAULT_SCOPE_PREFIX}${id}`;
    known.set(scope, id);
    return scope;
  };
  const threadOf =
    options.threadOf ??
    ((scope: string): ThreadId | undefined =>
      known.get(scope) ??
      (defaultScope && scope.startsWith(DEFAULT_SCOPE_PREFIX)
        ? threadId(scope.slice(DEFAULT_SCOPE_PREFIX.length))
        : undefined));
  const rows = Promise.resolve(options.rows);
  const loaded = new Map<string, ThreadState>();
  const listeners = new Set<(threadId: ThreadId) => void>();
  const remote = new Set<(threadId: ThreadId) => void>();
  let unwatch: Unsubscribe | undefined;
  let attaching = 0;

  const attach = () => {
    const attempt = ++attaching;
    void rows.then(
      (store) => {
        if (attempt !== attaching || remote.size === 0) return;
        unwatch = store.watch?.((scope) => {
          const id = threadOf(scope);
          if (id === undefined) return;
          for (const listener of [...remote]) {
            try {
              listener(id);
            } catch {
              // A listener cannot break the others.
            }
          }
        });
      },
      () => undefined,
    );
  };
  let tail: Promise<unknown> = Promise.resolve();

  const exclusive = <T>(work: () => Promise<T>): Promise<T> => {
    const run = tail.then(work, work);
    tail = run.catch(() => {});
    return run;
  };

  // A cached thread is used only while its rows are as it read them: another host on the same rows
  // (a second tab, another process) may have written since (cursor review F5, as the harness host).
  const load = async (id: ThreadId): Promise<ThreadState> => {
    const known = loaded.get(id);
    if (known && known.generation === (await (await rows).generation(scopeOf(id)))) return known;
    const state = stateOf(await (await rows).load(scopeOf(id)));
    loaded.set(id, state);
    return state;
  };

  /** A receipt for the command, in this thread or another one loaded. */
  const receiptIn = (draft: ThreadState, commandId: string): ThreadKernelReceipt | undefined => {
    const own = draft.receipts.get(commandId);
    if (own) return own;
    for (const state of loaded.values()) {
      const other = state.receipts.get(commandId);
      if (other) return other;
    }
    return undefined;
  };

  /**
   * Runs one command on a draft of its thread; when the command says its thread changed, writes the
   * difference through and only then keeps the draft.
   */
  const run = <T>(id: ThreadId, command: (draft: ThreadState, changed: () => void) => T): Promise<T> =>
    exclusive(async () => {
      // Written elsewhere since it was read: read again and run the command on that (bounded).
      for (let conflicts = 0; ; conflicts++) {
        const state = await load(id);
        const draft = draftOf(state);
        let changed = false;
        const value = command(draft, () => {
          changed = true;
        });
        if (!changed) return value;
        const writes = writesOf(state, draft);
        if (writes.length > 0) {
          try {
            await (await rows).apply(scopeOf(id), state.generation, writes);
          } catch (e) {
            if (e instanceof RowConflict && conflicts < MAX_CONFLICTS) {
              loaded.delete(id);
              continue;
            }
            throw e;
          }
          draft.generation = state.generation + 1;
        }
        loaded.set(id, draft);
        for (const listener of listeners) {
          try {
            listener(id);
          } catch {
            // A listener cannot turn a written command into a failure.
          }
        }
        return value;
      }
    });

  return {
    create: (input) =>
      run(input.threadId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const decision = prepareThreadCreate({
          command: input,
          ...defined('existingThread', s.thread),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          s.thread = {
            ...decision.result.value,
            limits: createThreadLimits(input.limits),
            activeStreamCount: 0,
          };
          remember(
            s,
            input.commandId,
            'create',
            { threadId: input.threadId, ownerId: input.ownerId, limits: input.limits },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    append: (input) =>
      run(input.threadId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const decision = prepareThreadAppendSettled({
          command: input,
          now,
          ...defined('thread', s.thread),
          ...defined('existingEntry', s.entries.get(input.entryId)),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          const { thread, entry } = decision.result.value;
          if (s.thread) s.thread = { ...s.thread, nextPosition: thread.nextPosition, lifecycle: thread.lifecycle };
          addEntry(s, entry);
          remember(
            s,
            input.commandId,
            'append',
            { threadId: input.threadId, entryId: input.entryId, descriptor: input.descriptor, payload: input.payload },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    begin: (input) =>
      run(input.threadId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const decision = prepareThreadBegin({
          command: input,
          now,
          ...defined('thread', s.thread),
          ...defined('existingEntry', s.entries.get(input.entryId)),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          const { thread, entry } = decision.result.value;
          if (s.thread)
            s.thread = {
              ...s.thread,
              nextPosition: thread.nextPosition,
              lifecycle: thread.lifecycle,
              activeStreamCount: s.thread.activeStreamCount + 1,
            };
          addEntry(s, entry);
          s.streams.set(input.entryId, {
            entryId: input.entryId,
            streamId: `thread:${input.threadId}:entry:${input.entryId}`,
            nextOffset: 0,
            lastActivityAt: now,
            retainedBytes: 0,
            retainedBatches: 0,
            status: 'active',
          });
          remember(
            s,
            input.commandId,
            'begin',
            { threadId: input.threadId, entryId: input.entryId, descriptor: input.descriptor },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    appendDelta: (input) =>
      run(input.threadId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const stored = s.streams.get(input.entryId);
        const decision = prepareThreadAppendDelta({
          command: input,
          now,
          ...defined('thread', s.thread),
          ...defined('entry', s.entries.get(input.entryId)),
          ...defined('stream', stored && kernelStream(stored)),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          const value = decision.result.value;
          s.entries.set(input.entryId, value.entry);
          if (value.kind === 'appended' && stored) {
            s.streams.set(input.entryId, {
              ...stored,
              nextOffset: value.nextOffset,
              lastActivityAt: now,
              retainedBytes: stored.retainedBytes + threadStoredEnvelopeBytes(input.payload),
              retainedBatches: stored.retainedBatches + 1,
            });
            s.deltas.set(stored.streamId, [
              ...(s.deltas.get(stored.streamId) ?? []),
              cloneData({
                entryId: input.entryId,
                start: input.start,
                endExclusive: value.nextOffset,
                payload: input.payload,
              }),
            ]);
          } else if (value.kind === 'aborted' && stored && s.thread) {
            s.thread = { ...s.thread, activeStreamCount: Math.max(0, s.thread.activeStreamCount - 1) };
            retire(s, s.thread, stored, now);
          }
          remember(
            s,
            input.commandId,
            'append-delta',
            { threadId: input.threadId, entryId: input.entryId, start: input.start, payload: input.payload },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    heartbeat: (input) =>
      run(input.threadId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const stored = s.streams.get(input.entryId);
        const decision = prepareThreadHeartbeat({
          command: input,
          now,
          ...defined('thread', s.thread),
          ...defined('entry', s.entries.get(input.entryId)),
          ...defined('stream', stored && kernelStream(stored)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && stored) {
          s.streams.set(input.entryId, { ...stored, lastActivityAt: Math.max(stored.lastActivityAt, now) });
          s.entries.set(input.entryId, decision.result.entry);
          changed();
        }
        return cloneData(decision.result);
      }),

    settle: (input) =>
      run(input.threadId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const stored = s.streams.get(input.entryId);
        const decision = prepareThreadSettle({
          command: input,
          now,
          ...defined('thread', s.thread),
          ...defined('entry', s.entries.get(input.entryId)),
          ...defined('stream', stored && kernelStream(stored)),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          if (s.thread) s.thread = { ...s.thread, activeStreamCount: Math.max(0, s.thread.activeStreamCount - 1) };
          s.entries.set(input.entryId, decision.result.value.entry);
          if (stored) {
            s.streams.delete(input.entryId);
            s.deltas.delete(stored.streamId);
          }
          remember(
            s,
            input.commandId,
            'settle',
            { threadId: input.threadId, entryId: input.entryId, payload: input.payload },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    abort: (input) =>
      run(input.threadId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const stored = s.streams.get(input.entryId);
        const decision = prepareThreadAbort({
          command: input,
          now,
          ...defined('thread', s.thread),
          ...defined('entry', s.entries.get(input.entryId)),
          ...defined('stream', stored && kernelStream(stored)),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          if (s.thread) {
            s.thread = { ...s.thread, activeStreamCount: Math.max(0, s.thread.activeStreamCount - 1) };
            if (stored) retire(s, s.thread, stored, now);
          }
          s.entries.set(input.entryId, decision.result.value.entry);
          remember(
            s,
            input.commandId,
            'abort',
            { threadId: input.threadId, entryId: input.entryId, reason: input.reason },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    window: (input) =>
      run(input.threadId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const thread = s.thread;
        if (!thread) return MISSING;
        const start = threadWindowRangeStart(input, thread.nextPosition);
        const take = threadPageTakeBound(input.limit, thread.limits.maxWindowEntries);
        const ids =
          input.direction === 'forward'
            ? s.order.slice(start, start + take)
            : s.order.slice(Math.max(0, start - take), start);
        const entries = ids.map((id) => s.entries.get(id)).filter((e): e is ThreadEntry => e !== undefined);
        const streams = new Map<string, ThreadKernelStream>();
        for (const e of entries) {
          const stored = s.streams.get(e.entryId);
          if (stored) streams.set(e.entryId, kernelStream(stored));
        }
        const projected = projectThreadWindow({ thread, entries, streams, request: input, now });
        return cloneData(projected);
      }),

    deltas: (input) =>
      run(input.threadId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const thread = s.thread;
        if (!thread) return MISSING;
        const stored = s.streams.get(input.entryId);
        const projected = projectThreadDeltas({
          thread,
          ...defined('entry', s.entries.get(input.entryId)),
          ...defined('stream', stored && kernelStream(stored)),
          deltas: deltaPage(
            s,
            stored,
            input.afterOffset,
            threadPageTakeBound(input.limit, thread.limits.maxDeltaPageBatches),
          ),
          request: input,
        });
        return cloneData(projected);
      }),

    remove: (input) =>
      run(input.threadId, (s, changed) => {
        if (applyInactivity(s, clock.now())) changed();
        const decision = prepareThreadRemove({
          command: input,
          ...defined('thread', s.thread),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          // The thread is fenced (removed) and everything else of it goes, its receipts too.
          s.entries.clear();
          s.order = [];
          s.streams.clear();
          s.deltas.clear();
          s.receipts.clear();
          if (s.thread) s.thread = { ...s.thread, lifecycle: 'removed', activeStreamCount: 0 };
          changed();
          const removed: ThreadCommandResult<ThreadRemoveResult> = {
            status: 'accepted',
            commandId: input.commandId,
            value: { threadId: input.threadId, lifecycle: 'removed' },
          };
          return cloneData(removed);
        }
        return cloneData(decision.result);
      }),

    subscribe(listener) {
      listeners.add(listener);
      return () => {
        listeners.delete(listener);
      };
    },

    watch(listener) {
      remote.add(listener);
      if (remote.size === 1) attach();
      return () => {
        if (!remote.delete(listener) || remote.size > 0) return;
        attaching++;
        unwatch?.();
        unwatch = undefined;
      };
    },

    forget: (id) =>
      exclusive(async () => {
        await (await rows).drop(scopeOf(id));
        loaded.delete(id);
      }),
  };
}

function defined<K extends string, T>(key: K, value: T | undefined): { readonly [P in K]: T } | Record<string, never> {
  return value === undefined ? {} : ({ [key]: value } as { readonly [P in K]: T });
}

function kernelStream(stream: Stream): ThreadKernelStream {
  const { entryId, nextOffset, lastActivityAt, retainedBytes, retainedBatches, status } = stream;
  return { entryId, nextOffset, lastActivityAt, retainedBytes, retainedBatches, status };
}

function addEntry(s: ThreadState, entry: ThreadEntry): void {
  s.entries.set(entry.entryId, entry);
  s.order.push(entry.entryId);
}

function remember(
  s: ThreadState,
  commandId: ThreadKernelReceipt['commandId'],
  kind: Parameters<typeof threadCommandFingerprint>[0],
  input: Readonly<Record<string, unknown>>,
  result: ThreadCommandResult<unknown>,
): void {
  s.receipts.set(commandId, {
    commandId,
    fingerprint: threadCommandFingerprint(kind, input),
    result: cloneData(result),
  });
}

/** A stream that ended without settling: kept a while as an aborted preview, or gone at once. */
function retire(s: ThreadState, thread: ThreadKernelThread, stored: Stream, now: number): void {
  if (thread.limits.abortedPreviewRetentionMs === 0) {
    s.streams.delete(stored.entryId);
    s.deltas.delete(stored.streamId);
    return;
  }
  s.streams.set(stored.entryId, {
    ...stored,
    status: 'aborted-preview',
    expiresAt: now + thread.limits.abortedPreviewRetentionMs,
  });
}

/** Drops expired previews and aborts streams idle past the timeout; true when anything changed. */
function applyInactivity(s: ThreadState, now: number): boolean {
  let changed = false;
  for (const stream of [...s.streams.values()]) {
    if (stream.status !== 'aborted-preview' || stream.expiresAt === undefined || stream.expiresAt > now) continue;
    s.streams.delete(stream.entryId);
    s.deltas.delete(stream.streamId);
    changed = true;
  }
  if (!s.thread) return changed;
  for (const stream of [...s.streams.values()]) {
    if (stream.status !== 'active') continue;
    const entry = s.entries.get(stream.entryId);
    if (!entry) continue;
    const decision = prepareThreadInactivityAbort({ now, thread: s.thread, entry, stream: kernelStream(stream) });
    if (decision.kind !== 'abort') continue;
    s.entries.set(stream.entryId, decision.entry);
    retire(s, s.thread, stream, now);
    s.thread = { ...s.thread, activeStreamCount: Math.max(0, s.thread.activeStreamCount - 1) };
    changed = true;
  }
  return changed;
}

/** The delta batches from `afterOffset`: the one covering it, then up to `take` after. */
function deltaPage(s: ThreadState, stored: Stream | undefined, afterOffset: number, take: number): ThreadDeltaBatch[] {
  if (!stored) return [];
  const all = s.deltas.get(stored.streamId) ?? [];
  let index = 0;
  while (index < all.length && (all[index] as ThreadDeltaBatch).start < afterOffset) index++;
  const covering = index > 0 ? all[index - 1] : undefined;
  const following = all.slice(index, index + take);
  return covering && covering.endExclusive > afterOffset ? [covering, ...following] : following;
}
