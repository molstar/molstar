/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { assertStrictJsonData, canonicalData } from './data.js';
import type {
  AbortThreadEntryCommand,
  AppendThreadDeltaCommand,
  AppendSettledThreadEntryCommand,
  BeginStreamingThreadEntryCommand,
  CreateThreadCommand,
  ThreadAppendDeltaResult,
  ThreadCommandResult,
  ThreadDeltaBatch,
  ThreadDeltaPage,
  ThreadDeltaResult,
  ThreadEntry,
  ThreadFailure,
  ThreadPayloadEnvelope,
  ThreadRemoveResult,
  ThreadSnapshot,
  ThreadWindow,
  ThreadWriteResult,
  HeartbeatThreadEntryCommand,
  ReadThreadDeltasRequest,
  ReadThreadWindowRequest,
  RemoveThreadCommand,
  SettleThreadEntryCommand,
} from './domain.js';
import { THREAD_PAYLOAD_VERSION } from './domain.js';
import type { ThreadCommandId, ThreadEntryId } from './ids.js';
import { createThreadLimits, type ThreadLimits } from './limits.js';

/**
 * The threads' kernel (analysis/viewer-mode.md §6): pure functions that decide a command against the state the host
 * gives them, and project windows and delta pages from it.
 */
export type ThreadCommandKind = 'create' | 'append' | 'begin' | 'append-delta' | 'settle' | 'abort' | 'remove';

export type ThreadKernelReceipt = Readonly<{
  readonly commandId: ThreadCommandId;
  readonly fingerprint: string;
  readonly result: ThreadCommandResult<unknown>;
}>;

export type ThreadKernelThread = ThreadSnapshot &
  Readonly<{
    readonly limits: ThreadLimits;
    readonly activeStreamCount: number;
  }>;

export type ThreadKernelStream = Readonly<{
  readonly entryId: ThreadEntryId;
  readonly nextOffset: number;
  readonly lastActivityAt: number;
  readonly retainedBytes: number;
  readonly retainedBatches: number;
  readonly status: 'active' | 'aborted-preview';
}>;

export type ThreadKernelDecision<T> =
  | Readonly<{ readonly kind: 'result'; readonly persist: boolean; readonly result: ThreadCommandResult<T> }>
  | Readonly<{ readonly kind: 'failure'; readonly failure: ThreadFailure }>;

export type ThreadKernelTouchDecision =
  | Readonly<{ readonly kind: 'result'; readonly persist: boolean; readonly result: ThreadDeltaResult }>
  | Readonly<{ readonly kind: 'failure'; readonly failure: ThreadFailure }>;

const CODEC_MAX_BYTES = 128;

function accepted<T>(commandId: ThreadCommandId, value: T): ThreadCommandResult<T> {
  return { status: 'accepted', commandId, value };
}

function failure(status: ThreadFailure['status'], reason: string): ThreadFailure {
  return { status, reason };
}

function replayOrConflict<T>(
  receipt: ThreadKernelReceipt | undefined,
  fingerprint: string,
): ThreadCommandResult<T> | ThreadFailure | undefined {
  if (receipt === undefined) return undefined;
  if (receipt.fingerprint !== fingerprint) return failure('conflict', 'command-reused');
  return receipt.result as ThreadCommandResult<T>;
}

export function threadStoredEnvelopeBytes(envelope: ThreadPayloadEnvelope): number {
  return new TextEncoder().encode(canonicalData(envelope)).byteLength;
}

export function validThreadPayloadEnvelope(value: unknown, maxBytes: number): value is ThreadPayloadEnvelope {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) return false;
  const envelope = value as Readonly<Record<string, unknown>>;
  if (Object.keys(envelope).length !== 3) return false;
  if (typeof envelope.codec !== 'string' || envelope.codec.trim().length === 0) return false;
  if (new TextEncoder().encode(envelope.codec).byteLength > CODEC_MAX_BYTES) return false;
  if (envelope.version !== THREAD_PAYLOAD_VERSION) return false;
  try {
    assertStrictJsonData(envelope.value);
  } catch {
    return false;
  }
  return threadStoredEnvelopeBytes(envelope as ThreadPayloadEnvelope) <= maxBytes;
}

export function threadCommandFingerprint(kind: ThreadCommandKind, input: Readonly<Record<string, unknown>>): string {
  const { commandId: _commandId, ...rest } = input;
  return canonicalData({ kind, ...rest });
}

function activeThread(thread: ThreadKernelThread | undefined): ThreadKernelThread | ThreadFailure {
  if (thread === undefined) return failure('not-found', 'thread-missing');
  if (thread.lifecycle !== 'active') return failure('closed', 'thread-closed');
  return thread;
}

function snapshot(thread: ThreadKernelThread): ThreadSnapshot {
  return {
    threadId: thread.threadId,
    ownerId: thread.ownerId,
    nextPosition: thread.nextPosition,
    lifecycle: thread.lifecycle,
  };
}

function projectStreaming(entry: ThreadEntry, stream: ThreadKernelStream, now: number): ThreadEntry {
  return {
    ...entry,
    state: { kind: 'streaming', nextOffset: stream.nextOffset, lastActivityAt: stream.lastActivityAt ?? now },
  };
}

export function prepareThreadCreate(input: {
  readonly command: CreateThreadCommand;
  readonly existingThread?: ThreadKernelThread;
  readonly existingReceipt?: ThreadKernelReceipt;
}): ThreadKernelDecision<ThreadSnapshot> {
  const fingerprint = threadCommandFingerprint('create', {
    threadId: input.command.threadId,
    ownerId: input.command.ownerId,
    limits: input.command.limits,
  });
  const replay = replayOrConflict<ThreadSnapshot>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as ThreadFailure };
  }
  if (input.existingThread !== undefined) return { kind: 'failure', failure: failure('conflict', 'thread-exists') };
  let limits: ThreadLimits;
  try {
    limits = createThreadLimits(input.command.limits);
  } catch {
    return { kind: 'failure', failure: failure('invalid', 'limits') };
  }
  const thread: ThreadKernelThread = {
    threadId: input.command.threadId,
    ownerId: input.command.ownerId,
    nextPosition: 0,
    lifecycle: 'active',
    limits,
    activeStreamCount: 0,
  };
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, snapshot(thread)),
  };
}

export function prepareThreadAppendSettled(input: {
  readonly command: AppendSettledThreadEntryCommand;
  readonly now: number;
  readonly thread?: ThreadKernelThread;
  readonly existingEntry?: ThreadEntry;
  readonly existingReceipt?: ThreadKernelReceipt;
}): ThreadKernelDecision<ThreadWriteResult> {
  const fingerprint = threadCommandFingerprint('append', {
    threadId: input.command.threadId,
    entryId: input.command.entryId,
    descriptor: input.command.descriptor,
    payload: input.command.payload,
  });
  const replay = replayOrConflict<ThreadWriteResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as ThreadFailure };
  }
  const thread = activeThread(input.thread);
  if ('status' in thread) return { kind: 'failure', failure: thread };
  if (!validThreadPayloadEnvelope(input.command.descriptor, thread.limits.maxEntryBytes)) {
    return { kind: 'failure', failure: failure('invalid', 'descriptor') };
  }
  if (!validThreadPayloadEnvelope(input.command.payload, thread.limits.maxEntryBytes)) {
    return { kind: 'failure', failure: failure('invalid', 'payload') };
  }
  if (input.existingEntry !== undefined) return { kind: 'failure', failure: failure('conflict', 'entry-exists') };
  const entry: ThreadEntry = {
    entryId: input.command.entryId,
    threadId: thread.threadId,
    position: thread.nextPosition,
    descriptor: input.command.descriptor,
    state: { kind: 'settled', payload: input.command.payload, settledAt: input.now },
  };
  const next: ThreadKernelThread = { ...thread, nextPosition: thread.nextPosition + 1 };
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, { thread: snapshot(next), entry }),
  };
}

export function prepareThreadBegin(input: {
  readonly command: BeginStreamingThreadEntryCommand;
  readonly now: number;
  readonly thread?: ThreadKernelThread;
  readonly existingEntry?: ThreadEntry;
  readonly existingReceipt?: ThreadKernelReceipt;
}): ThreadKernelDecision<ThreadWriteResult> {
  const fingerprint = threadCommandFingerprint('begin', {
    threadId: input.command.threadId,
    entryId: input.command.entryId,
    descriptor: input.command.descriptor,
  });
  const replay = replayOrConflict<ThreadWriteResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as ThreadFailure };
  }
  const thread = activeThread(input.thread);
  if ('status' in thread) return { kind: 'failure', failure: thread };
  if (thread.activeStreamCount >= thread.limits.maxActiveStreams) {
    return { kind: 'failure', failure: failure('invalid', 'active-stream-limit') };
  }
  if (!validThreadPayloadEnvelope(input.command.descriptor, thread.limits.maxEntryBytes)) {
    return { kind: 'failure', failure: failure('invalid', 'descriptor') };
  }
  if (input.existingEntry !== undefined) return { kind: 'failure', failure: failure('conflict', 'entry-exists') };
  const stream: ThreadKernelStream = {
    entryId: input.command.entryId,
    nextOffset: 0,
    lastActivityAt: input.now,
    retainedBytes: 0,
    retainedBatches: 0,
    status: 'active',
  };
  const entry: ThreadEntry = {
    entryId: input.command.entryId,
    threadId: thread.threadId,
    position: thread.nextPosition,
    descriptor: input.command.descriptor,
    state: { kind: 'streaming', nextOffset: 0, lastActivityAt: input.now },
  };
  const next: ThreadKernelThread = {
    ...thread,
    nextPosition: thread.nextPosition + 1,
    activeStreamCount: thread.activeStreamCount + 1,
  };
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, {
      thread: snapshot(next),
      entry: projectStreaming(entry, stream, input.now),
    }),
  };
}

export function prepareThreadAppendDelta(input: {
  readonly command: AppendThreadDeltaCommand;
  readonly now: number;
  readonly thread?: ThreadKernelThread;
  readonly entry?: ThreadEntry;
  readonly stream?: ThreadKernelStream;
  readonly existingReceipt?: ThreadKernelReceipt;
}): ThreadKernelDecision<ThreadAppendDeltaResult> {
  const fingerprint = threadCommandFingerprint('append-delta', {
    threadId: input.command.threadId,
    entryId: input.command.entryId,
    start: input.command.start,
    payload: input.command.payload,
  });
  const replay = replayOrConflict<ThreadAppendDeltaResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as ThreadFailure };
  }
  const thread = activeThread(input.thread);
  if ('status' in thread) return { kind: 'failure', failure: thread };
  if (input.entry === undefined || input.stream === undefined || input.stream.status !== 'active') {
    return { kind: 'failure', failure: failure(input.entry === undefined ? 'not-found' : 'invalid', 'not-streaming') };
  }
  if (input.entry.state.kind !== 'streaming') {
    return { kind: 'failure', failure: failure('invalid', 'terminal-entry') };
  }
  if (!validThreadPayloadEnvelope(input.command.payload, thread.limits.maxBatchBytes)) {
    return { kind: 'failure', failure: failure('invalid', 'payload') };
  }
  if (input.command.start !== input.stream.nextOffset) {
    return { kind: 'failure', failure: failure('invalid', 'delta-offset') };
  }
  const batchBytes = threadStoredEnvelopeBytes(input.command.payload);
  if (
    input.stream.retainedBytes + batchBytes > thread.limits.maxActiveStreamBytes ||
    input.stream.retainedBatches + 1 > thread.limits.maxActiveStreamBatches
  ) {
    const entry: ThreadEntry = {
      ...input.entry,
      state: { kind: 'aborted', reason: 'limit-exceeded', endedAt: input.now },
    };
    const next: ThreadKernelThread = { ...thread, activeStreamCount: Math.max(0, thread.activeStreamCount - 1) };
    return {
      kind: 'result',
      persist: true,
      result: accepted(input.command.commandId, {
        kind: 'aborted',
        thread: snapshot(next),
        entry,
        reason: 'limit-exceeded',
      }),
    };
  }
  const endExclusive = input.command.start + batchBytes;
  const stream: ThreadKernelStream = {
    ...input.stream,
    nextOffset: endExclusive,
    lastActivityAt: input.now,
    retainedBytes: input.stream.retainedBytes + batchBytes,
    retainedBatches: input.stream.retainedBatches + 1,
  };
  const entry = projectStreaming(input.entry, stream, input.now);
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, {
      kind: 'appended',
      thread: snapshot(thread),
      entry,
      nextOffset: stream.nextOffset,
    }),
  };
}

export function prepareThreadHeartbeat(input: {
  readonly command: HeartbeatThreadEntryCommand;
  readonly now: number;
  readonly thread?: ThreadKernelThread;
  readonly entry?: ThreadEntry;
  readonly stream?: ThreadKernelStream;
}): ThreadKernelTouchDecision {
  const thread = activeThread(input.thread);
  if ('status' in thread) return { kind: 'failure', failure: thread };
  if (input.entry === undefined || input.stream === undefined || input.stream.status !== 'active') {
    return { kind: 'failure', failure: failure(input.entry === undefined ? 'not-found' : 'invalid', 'not-streaming') };
  }
  const lastActivityAt = Math.max(input.stream.lastActivityAt, input.now);
  const stream: ThreadKernelStream = { ...input.stream, lastActivityAt };
  const entry = projectStreaming(input.entry, stream, input.now);
  return {
    kind: 'result',
    persist: lastActivityAt > input.stream.lastActivityAt,
    result: {
      thread: snapshot(thread),
      entry,
      nextOffset: stream.nextOffset,
    },
  };
}

export function prepareThreadSettle(input: {
  readonly command: SettleThreadEntryCommand;
  readonly now: number;
  readonly thread?: ThreadKernelThread;
  readonly entry?: ThreadEntry;
  readonly stream?: ThreadKernelStream;
  readonly existingReceipt?: ThreadKernelReceipt;
}): ThreadKernelDecision<ThreadWriteResult> {
  const fingerprint = threadCommandFingerprint('settle', {
    threadId: input.command.threadId,
    entryId: input.command.entryId,
    payload: input.command.payload,
  });
  const replay = replayOrConflict<ThreadWriteResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as ThreadFailure };
  }
  const thread = activeThread(input.thread);
  if ('status' in thread) return { kind: 'failure', failure: thread };
  if (input.entry === undefined) return { kind: 'failure', failure: failure('not-found', 'entry-missing') };
  if (input.entry.state.kind === 'settled') {
    return { kind: 'failure', failure: failure('invalid', 'already-settled') };
  }
  if (input.entry.state.kind === 'aborted' || input.stream === undefined || input.stream.status !== 'active') {
    return { kind: 'failure', failure: failure('invalid', 'terminal-entry') };
  }
  if (!validThreadPayloadEnvelope(input.command.payload, thread.limits.maxEntryBytes)) {
    return { kind: 'failure', failure: failure('invalid', 'payload') };
  }
  const entry: ThreadEntry = {
    ...input.entry,
    state: { kind: 'settled', payload: input.command.payload, settledAt: input.now },
  };
  const next: ThreadKernelThread = { ...thread, activeStreamCount: thread.activeStreamCount - 1 };
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, { thread: snapshot(next), entry }),
  };
}

export function prepareThreadAbort(input: {
  readonly command: AbortThreadEntryCommand;
  readonly now: number;
  readonly thread?: ThreadKernelThread;
  readonly entry?: ThreadEntry;
  readonly stream?: ThreadKernelStream;
  readonly existingReceipt?: ThreadKernelReceipt;
}): ThreadKernelDecision<ThreadWriteResult> {
  const fingerprint = threadCommandFingerprint('abort', {
    threadId: input.command.threadId,
    entryId: input.command.entryId,
    reason: input.command.reason,
  });
  const replay = replayOrConflict<ThreadWriteResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as ThreadFailure };
  }
  const thread = activeThread(input.thread);
  if ('status' in thread) return { kind: 'failure', failure: thread };
  if (input.entry === undefined) return { kind: 'failure', failure: failure('not-found', 'entry-missing') };
  if (input.entry.state.kind !== 'streaming' || input.stream === undefined || input.stream.status !== 'active') {
    return { kind: 'failure', failure: failure('invalid', 'terminal-entry') };
  }
  const entry: ThreadEntry = {
    ...input.entry,
    state: { kind: 'aborted', reason: input.command.reason, endedAt: input.now },
  };
  const next: ThreadKernelThread = { ...thread, activeStreamCount: thread.activeStreamCount - 1 };
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, { thread: snapshot(next), entry }),
  };
}

export function prepareThreadRemove(input: {
  readonly command: RemoveThreadCommand;
  readonly thread?: ThreadKernelThread;
  readonly existingReceipt?: ThreadKernelReceipt;
}): ThreadKernelDecision<ThreadRemoveResult> {
  const fingerprint = threadCommandFingerprint('remove', { threadId: input.command.threadId });
  const replay = replayOrConflict<ThreadRemoveResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as ThreadFailure };
  }
  if (input.thread === undefined) {
    return {
      kind: 'result',
      persist: false,
      result: accepted(input.command.commandId, { threadId: input.command.threadId, lifecycle: 'removed' }),
    };
  }
  switch (input.thread.lifecycle) {
    case 'removed':
      return {
        kind: 'result',
        persist: false,
        result: accepted(input.command.commandId, { threadId: input.command.threadId, lifecycle: 'removed' }),
      };
    case 'closing':
      return {
        kind: 'result',
        persist: false,
        result: accepted(input.command.commandId, { threadId: input.command.threadId, lifecycle: 'closing' }),
      };
    case 'active':
      return {
        kind: 'result',
        persist: true,
        result: accepted(input.command.commandId, { threadId: input.command.threadId, lifecycle: 'closing' }),
      };
    default: {
      const exhaustive: never = input.thread.lifecycle;
      throw new Error(`unexpected thread lifecycle: ${exhaustive}`);
    }
  }
}

export function prepareThreadInactivityAbort(input: {
  readonly now: number;
  readonly thread: ThreadKernelThread;
  readonly entry: ThreadEntry;
  readonly stream: ThreadKernelStream;
}):
  | Readonly<{ readonly kind: 'skip' }>
  | Readonly<{ readonly kind: 'abort'; readonly thread: ThreadSnapshot; readonly entry: ThreadEntry }> {
  if (input.stream.status !== 'active' || input.entry.state.kind !== 'streaming') return { kind: 'skip' };
  if (input.now - input.stream.lastActivityAt < input.thread.limits.inactivityTimeoutMs) return { kind: 'skip' };
  return {
    kind: 'abort',
    thread: snapshot({ ...input.thread, activeStreamCount: Math.max(0, input.thread.activeStreamCount - 1) }),
    entry: {
      ...input.entry,
      state: { kind: 'aborted', reason: 'inactivity', endedAt: input.now },
    },
  };
}

function entryWindowBytes(entry: ThreadEntry): number {
  const descriptor = threadStoredEnvelopeBytes(entry.descriptor);
  return entry.state.kind === 'settled' ? descriptor + threadStoredEnvelopeBytes(entry.state.payload) : descriptor;
}

export function threadPageTakeBound(limit: number, ceiling: number): number {
  return Math.min(Math.max(limit, 0), ceiling) + 1;
}

export function threadWindowRangeStart(
  request: Pick<ReadThreadWindowRequest, 'start' | 'direction'>,
  nextPosition: number,
): number {
  if (request.start !== undefined) return request.start;
  return request.direction === 'forward' ? 0 : nextPosition;
}

export function projectThreadWindow(input: {
  readonly thread: ThreadKernelThread;
  readonly entries: readonly ThreadEntry[];
  readonly streams: ReadonlyMap<string, ThreadKernelStream>;
  readonly request: ReadThreadWindowRequest;
  readonly now: number;
}): ThreadWindow | ThreadFailure {
  if (input.request.limit < 1 || input.request.limit > input.thread.limits.maxWindowEntries) {
    return failure('invalid', 'window-limit');
  }
  const ordered = [...input.entries].sort((left, right) => left.position - right.position);
  const start = threadWindowRangeStart(input.request, input.thread.nextPosition);
  const forward = input.request.direction === 'forward';
  const candidates = forward
    ? ordered.filter((entry) => entry.position >= start)
    : [...ordered.filter((entry) => entry.position < start)].reverse();
  const projected: ThreadEntry[] = [];
  let bytes = 0;
  for (const entry of candidates) {
    const stream = input.streams.get(entry.entryId);
    const next =
      entry.state.kind === 'streaming' && stream !== undefined ? projectStreaming(entry, stream, input.now) : entry;
    const size = entryWindowBytes(next);
    if (projected.length >= input.request.limit) break;
    if (projected.length > 0 && bytes + size > input.thread.limits.maxWindowBytes) break;
    if (projected.length === 0 && size > input.thread.limits.maxWindowBytes) {
      return failure('invalid', 'window-bytes');
    }
    projected.push(next);
    bytes += size;
  }
  if (!forward) projected.reverse();
  const first = projected[0]?.position ?? start;
  const last = projected[projected.length - 1]?.position;
  return {
    threadId: input.thread.threadId,
    entries: projected,
    start: first,
    endExclusive: last === undefined ? first : last + 1,
    hasMore: candidates.length > projected.length,
  };
}

export function projectThreadDeltas(input: {
  readonly thread: ThreadKernelThread;
  readonly entry?: ThreadEntry;
  readonly stream?: ThreadKernelStream;
  readonly deltas: readonly ThreadDeltaBatch[];
  readonly request: ReadThreadDeltasRequest;
}): ThreadDeltaPage | ThreadFailure {
  if (input.entry === undefined) return failure('not-found', 'entry-missing');
  if (input.request.limit < 1 || input.request.limit > input.thread.limits.maxDeltaPageBatches) {
    return failure('invalid', 'delta-page-limit');
  }
  if (input.entry.state.kind !== 'streaming' || input.stream === undefined) {
    return {
      entryId: input.entry.entryId,
      batches: [],
      nextOffset: input.request.afterOffset,
      hasMore: false,
    };
  }
  const ordered = [...input.deltas].sort((left, right) => left.start - right.start);
  const next = ordered.find((batch) => batch.start >= input.request.afterOffset);
  if (input.request.afterOffset > 0 && (next === undefined || next.start !== input.request.afterOffset)) {
    const covering = ordered.find(
      (batch) => batch.start <= input.request.afterOffset && batch.endExclusive > input.request.afterOffset,
    );
    if (covering === undefined && next !== undefined && next.start > input.request.afterOffset) {
      return {
        entryId: input.entry.entryId,
        batches: [],
        nextOffset: input.request.afterOffset,
        hasMore: false,
        gap: { afterOffset: input.request.afterOffset },
      };
    }
  }
  const page: ThreadDeltaBatch[] = [];
  let offset = input.request.afterOffset;
  let bytes = 0;
  for (const batch of ordered) {
    if (batch.start < input.request.afterOffset) continue;
    if (batch.start !== offset) {
      if (page.length === 0) {
        return {
          entryId: input.entry.entryId,
          batches: [],
          nextOffset: input.request.afterOffset,
          hasMore: false,
          gap: { afterOffset: input.request.afterOffset },
        };
      }
      break;
    }
    const size = threadStoredEnvelopeBytes(batch.payload);
    if (page.length >= input.request.limit || bytes + size > input.thread.limits.maxDeltaPageBytes) break;
    page.push(batch);
    offset = batch.endExclusive;
    bytes += size;
  }
  return {
    entryId: input.entry.entryId,
    batches: page,
    nextOffset: offset,
    hasMore: ordered.some((batch) => batch.start >= offset),
  };
}
