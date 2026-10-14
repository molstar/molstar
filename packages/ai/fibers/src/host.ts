/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { cloneData } from './data.js';
import type {
  FiberClock,
  FiberCommandResult,
  FiberDeltaBatch,
  FiberEntry,
  FiberFailure,
  FiberPort,
  FiberRemoveResult,
} from './domain.js';
import { type FiberId, fiberId } from './ids.js';
import {
  type FiberKernelFiber,
  type FiberKernelReceipt,
  type FiberKernelStream,
  fiberCommandFingerprint,
  fiberPageTakeBound,
  fiberStoredEnvelopeBytes,
  fiberWindowRangeStart,
  prepareFiberAbort,
  prepareFiberAppendDelta,
  prepareFiberAppendSettled,
  prepareFiberBegin,
  prepareFiberCreate,
  prepareFiberHeartbeat,
  prepareFiberInactivityAbort,
  prepareFiberRemove,
  prepareFiberSettle,
  projectFiberDeltas,
  projectFiberWindow,
} from './kernel.js';
import { createFiberLimits } from './limits.js';
import { draftOf, type FiberState, type Stream, stateOf, writesOf } from './state.js';
import { type FiberRows, RowConflict, type Unsubscribe } from './storage.js';

/**
 * The fibers' host (analysis/viewer-mode.md §6), written from anyfold's in-memory host and its Convex host as
 * Stories ported it: the kernel decides every command; the host keeps each fiber's state and writes what a command
 * changed through to its rows before the command returns (so an accepted command is durable). Unlike the in-memory
 * host, a fiber is its own scope of rows, loaded when first used and kept in memory after.
 *
 * Commands run one at a time. A replayed command id is found in its fiber's receipts, and in the receipts of the
 * other fibers loaded here (a command id reused for another fiber is a conflict); callers name commands after
 * their fiber, so ids never meet across fibers.
 */
export interface FiberHost extends FiberPort {
  /** Hears which fiber a command of this host changed. */
  subscribe(listener: (fiberId: FiberId) => void): Unsubscribe;
  /**
   * Hears which fiber was changed through another handle to the same rows (another tab, another process), when
   * the rows can tell (`RowStore.watch`); otherwise never. The next command or read on that fiber sees the change
   * with or without this.
   */
  watch(listener: (fiberId: FiberId) => void): Unsubscribe;
  /** Forgets a fiber entirely: its rows and what is loaded (its owner was removed). */
  forget(fiberId: FiberId): Promise<void>;
}

export interface FiberHostOptions {
  rows: FiberRows | Promise<FiberRows>;
  clock?: FiberClock;
  /** A fiber's scope in the rows; default `fiber:<id>`. */
  scopeOf?: (fiberId: FiberId) => string;
  /**
   * The fiber of a scope the rows report as changed (`watch`); default: the fibers this host has used, and with
   * the default `scopeOf` any `fiber:<id>`.
   */
  fiberOf?: (scope: string) => FiberId | undefined;
}

const MAX_CONFLICTS = 8;
/** A read of a fiber this store does not hold. */
const MISSING: FiberFailure = { status: 'not-found', reason: 'fiber-missing' };
const DEFAULT_SCOPE_PREFIX = 'fiber:';

export function createFiberHost(options: FiberHostOptions): FiberHost {
  const clock = options.clock ?? { now: () => Date.now() };
  const defaultScope = options.scopeOf === undefined;
  const known = new Map<string, FiberId>();
  const scopeOf = (id: FiberId): string => {
    const scope = options.scopeOf?.(id) ?? `${DEFAULT_SCOPE_PREFIX}${id}`;
    known.set(scope, id);
    return scope;
  };
  const fiberOf =
    options.fiberOf ??
    ((scope: string): FiberId | undefined =>
      known.get(scope) ??
      (defaultScope && scope.startsWith(DEFAULT_SCOPE_PREFIX)
        ? fiberId(scope.slice(DEFAULT_SCOPE_PREFIX.length))
        : undefined));
  const rows = Promise.resolve(options.rows);
  const loaded = new Map<string, FiberState>();
  const listeners = new Set<(fiberId: FiberId) => void>();
  const remote = new Set<(fiberId: FiberId) => void>();
  let unwatch: Unsubscribe | undefined;
  let attaching = 0;

  const attach = () => {
    const attempt = ++attaching;
    void rows.then(
      (store) => {
        if (attempt !== attaching || remote.size === 0) return;
        unwatch = store.watch?.((scope) => {
          const id = fiberOf(scope);
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

  // A cached fiber is used only while its rows are as it read them: another host on the same rows
  // (a second tab, another process) may have written since (cursor review F5, as the harness host).
  const load = async (id: FiberId): Promise<FiberState> => {
    const known = loaded.get(id);
    if (known && known.generation === (await (await rows).generation(scopeOf(id)))) return known;
    const state = stateOf(await (await rows).load(scopeOf(id)));
    loaded.set(id, state);
    return state;
  };

  /** A receipt for the command, in this fiber or another one loaded. */
  const receiptIn = (draft: FiberState, commandId: string): FiberKernelReceipt | undefined => {
    const own = draft.receipts.get(commandId);
    if (own) return own;
    for (const state of loaded.values()) {
      const other = state.receipts.get(commandId);
      if (other) return other;
    }
    return undefined;
  };

  /**
   * Runs one command on a draft of its fiber; when the command says its fiber changed, writes the
   * difference through and only then keeps the draft.
   */
  const run = <T>(id: FiberId, command: (draft: FiberState, changed: () => void) => T): Promise<T> =>
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
      run(input.fiberId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const decision = prepareFiberCreate({
          command: input,
          ...defined('existingFiber', s.fiber),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          s.fiber = {
            ...decision.result.value,
            limits: createFiberLimits(input.limits),
            activeStreamCount: 0,
          };
          remember(
            s,
            input.commandId,
            'create',
            { fiberId: input.fiberId, ownerId: input.ownerId, limits: input.limits },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    append: (input) =>
      run(input.fiberId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const decision = prepareFiberAppendSettled({
          command: input,
          now,
          ...defined('fiber', s.fiber),
          ...defined('existingEntry', s.entries.get(input.entryId)),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          const { fiber, entry } = decision.result.value;
          if (s.fiber) s.fiber = { ...s.fiber, nextPosition: fiber.nextPosition, lifecycle: fiber.lifecycle };
          addEntry(s, entry);
          remember(
            s,
            input.commandId,
            'append',
            { fiberId: input.fiberId, entryId: input.entryId, descriptor: input.descriptor, payload: input.payload },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    begin: (input) =>
      run(input.fiberId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const decision = prepareFiberBegin({
          command: input,
          now,
          ...defined('fiber', s.fiber),
          ...defined('existingEntry', s.entries.get(input.entryId)),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          const { fiber, entry } = decision.result.value;
          if (s.fiber)
            s.fiber = {
              ...s.fiber,
              nextPosition: fiber.nextPosition,
              lifecycle: fiber.lifecycle,
              activeStreamCount: s.fiber.activeStreamCount + 1,
            };
          addEntry(s, entry);
          s.streams.set(input.entryId, {
            entryId: input.entryId,
            streamId: `fiber:${input.fiberId}:entry:${input.entryId}`,
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
            { fiberId: input.fiberId, entryId: input.entryId, descriptor: input.descriptor },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    appendDelta: (input) =>
      run(input.fiberId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const stored = s.streams.get(input.entryId);
        const decision = prepareFiberAppendDelta({
          command: input,
          now,
          ...defined('fiber', s.fiber),
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
              retainedBytes: stored.retainedBytes + fiberStoredEnvelopeBytes(input.payload),
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
          } else if (value.kind === 'aborted' && stored && s.fiber) {
            s.fiber = { ...s.fiber, activeStreamCount: Math.max(0, s.fiber.activeStreamCount - 1) };
            retire(s, s.fiber, stored, now);
          }
          remember(
            s,
            input.commandId,
            'append-delta',
            { fiberId: input.fiberId, entryId: input.entryId, start: input.start, payload: input.payload },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    heartbeat: (input) =>
      run(input.fiberId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const stored = s.streams.get(input.entryId);
        const decision = prepareFiberHeartbeat({
          command: input,
          now,
          ...defined('fiber', s.fiber),
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
      run(input.fiberId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const stored = s.streams.get(input.entryId);
        const decision = prepareFiberSettle({
          command: input,
          now,
          ...defined('fiber', s.fiber),
          ...defined('entry', s.entries.get(input.entryId)),
          ...defined('stream', stored && kernelStream(stored)),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          if (s.fiber) s.fiber = { ...s.fiber, activeStreamCount: Math.max(0, s.fiber.activeStreamCount - 1) };
          s.entries.set(input.entryId, decision.result.value.entry);
          if (stored) {
            s.streams.delete(input.entryId);
            s.deltas.delete(stored.streamId);
          }
          remember(
            s,
            input.commandId,
            'settle',
            { fiberId: input.fiberId, entryId: input.entryId, payload: input.payload },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    abort: (input) =>
      run(input.fiberId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const stored = s.streams.get(input.entryId);
        const decision = prepareFiberAbort({
          command: input,
          now,
          ...defined('fiber', s.fiber),
          ...defined('entry', s.entries.get(input.entryId)),
          ...defined('stream', stored && kernelStream(stored)),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          if (s.fiber) {
            s.fiber = { ...s.fiber, activeStreamCount: Math.max(0, s.fiber.activeStreamCount - 1) };
            if (stored) retire(s, s.fiber, stored, now);
          }
          s.entries.set(input.entryId, decision.result.value.entry);
          remember(
            s,
            input.commandId,
            'abort',
            { fiberId: input.fiberId, entryId: input.entryId, reason: input.reason },
            decision.result,
          );
          changed();
        }
        return cloneData(decision.result);
      }),

    window: (input) =>
      run(input.fiberId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const fiber = s.fiber;
        if (!fiber) return MISSING;
        const start = fiberWindowRangeStart(input, fiber.nextPosition);
        const take = fiberPageTakeBound(input.limit, fiber.limits.maxWindowEntries);
        const ids =
          input.direction === 'forward'
            ? s.order.slice(start, start + take)
            : s.order.slice(Math.max(0, start - take), start);
        const entries = ids.map((id) => s.entries.get(id)).filter((e): e is FiberEntry => e !== undefined);
        const streams = new Map<string, FiberKernelStream>();
        for (const e of entries) {
          const stored = s.streams.get(e.entryId);
          if (stored) streams.set(e.entryId, kernelStream(stored));
        }
        const projected = projectFiberWindow({ fiber, entries, streams, request: input, now });
        return cloneData(projected);
      }),

    deltas: (input) =>
      run(input.fiberId, (s, changed) => {
        const now = clock.now();
        if (applyInactivity(s, now)) changed();
        const fiber = s.fiber;
        if (!fiber) return MISSING;
        const stored = s.streams.get(input.entryId);
        const projected = projectFiberDeltas({
          fiber,
          ...defined('entry', s.entries.get(input.entryId)),
          ...defined('stream', stored && kernelStream(stored)),
          deltas: deltaPage(
            s,
            stored,
            input.afterOffset,
            fiberPageTakeBound(input.limit, fiber.limits.maxDeltaPageBatches),
          ),
          request: input,
        });
        return cloneData(projected);
      }),

    remove: (input) =>
      run(input.fiberId, (s, changed) => {
        if (applyInactivity(s, clock.now())) changed();
        const decision = prepareFiberRemove({
          command: input,
          ...defined('fiber', s.fiber),
          ...defined('existingReceipt', receiptIn(s, input.commandId)),
        });
        if (decision.kind === 'failure') return decision.failure;
        if (decision.persist && decision.result.status === 'accepted') {
          // The fiber is fenced (removed) and everything else of it goes, its receipts too.
          s.entries.clear();
          s.order = [];
          s.streams.clear();
          s.deltas.clear();
          s.receipts.clear();
          if (s.fiber) s.fiber = { ...s.fiber, lifecycle: 'removed', activeStreamCount: 0 };
          changed();
          const removed: FiberCommandResult<FiberRemoveResult> = {
            status: 'accepted',
            commandId: input.commandId,
            value: { fiberId: input.fiberId, lifecycle: 'removed' },
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

function kernelStream(stream: Stream): FiberKernelStream {
  const { entryId, nextOffset, lastActivityAt, retainedBytes, retainedBatches, status } = stream;
  return { entryId, nextOffset, lastActivityAt, retainedBytes, retainedBatches, status };
}

function addEntry(s: FiberState, entry: FiberEntry): void {
  s.entries.set(entry.entryId, entry);
  s.order.push(entry.entryId);
}

function remember(
  s: FiberState,
  commandId: FiberKernelReceipt['commandId'],
  kind: Parameters<typeof fiberCommandFingerprint>[0],
  input: Readonly<Record<string, unknown>>,
  result: FiberCommandResult<unknown>,
): void {
  s.receipts.set(commandId, {
    commandId,
    fingerprint: fiberCommandFingerprint(kind, input),
    result: cloneData(result),
  });
}

/** A stream that ended without settling: kept a while as an aborted preview, or gone at once. */
function retire(s: FiberState, fiber: FiberKernelFiber, stored: Stream, now: number): void {
  if (fiber.limits.abortedPreviewRetentionMs === 0) {
    s.streams.delete(stored.entryId);
    s.deltas.delete(stored.streamId);
    return;
  }
  s.streams.set(stored.entryId, {
    ...stored,
    status: 'aborted-preview',
    expiresAt: now + fiber.limits.abortedPreviewRetentionMs,
  });
}

/** Drops expired previews and aborts streams idle past the timeout; true when anything changed. */
function applyInactivity(s: FiberState, now: number): boolean {
  let changed = false;
  for (const stream of [...s.streams.values()]) {
    if (stream.status !== 'aborted-preview' || stream.expiresAt === undefined || stream.expiresAt > now) continue;
    s.streams.delete(stream.entryId);
    s.deltas.delete(stream.streamId);
    changed = true;
  }
  if (!s.fiber) return changed;
  for (const stream of [...s.streams.values()]) {
    if (stream.status !== 'active') continue;
    const entry = s.entries.get(stream.entryId);
    if (!entry) continue;
    const decision = prepareFiberInactivityAbort({ now, fiber: s.fiber, entry, stream: kernelStream(stream) });
    if (decision.kind !== 'abort') continue;
    s.entries.set(stream.entryId, decision.entry);
    retire(s, s.fiber, stream, now);
    s.fiber = { ...s.fiber, activeStreamCount: Math.max(0, s.fiber.activeStreamCount - 1) };
    changed = true;
  }
  return changed;
}

/** The delta batches from `afterOffset`: the one covering it, then up to `take` after. */
function deltaPage(s: FiberState, stored: Stream | undefined, afterOffset: number, take: number): FiberDeltaBatch[] {
  if (!stored) return [];
  const all = s.deltas.get(stored.streamId) ?? [];
  let index = 0;
  while (index < all.length && (all[index] as FiberDeltaBatch).start < afterOffset) index++;
  const covering = index > 0 ? all[index - 1] : undefined;
  const following = all.slice(index, index + take);
  return covering && covering.endExclusive > afterOffset ? [covering, ...following] : following;
}
