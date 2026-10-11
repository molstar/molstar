/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { assertStrictJsonData, canonicalData } from './data.js';
import type {
  AbortFiberEntryCommand,
  AppendFiberDeltaCommand,
  AppendSettledFiberEntryCommand,
  BeginStreamingFiberEntryCommand,
  CreateFiberCommand,
  FiberAppendDeltaResult,
  FiberCleanupCursor,
  FiberCleanupPage,
  FiberCommandResult,
  FiberDeltaBatch,
  FiberDeltaPage,
  FiberDeltaResult,
  FiberEntry,
  FiberFailure,
  FiberPayloadEnvelope,
  FiberRemoveResult,
  FiberSnapshot,
  FiberWindow,
  FiberWriteResult,
  HeartbeatFiberEntryCommand,
  ReadFiberDeltasRequest,
  ReadFiberWindowRequest,
  RemoveFiberCommand,
  SettleFiberEntryCommand,
} from './domain.js';
import { FIBER_PAYLOAD_VERSION } from './domain.js';
import type { FiberCommandId, FiberEntryId } from './ids.js';
import { createFiberLimits, type FiberLimits } from './limits.js';

/**
 * The fibers' kernel (analysis/viewer-mode.md §6): pure functions that decide a command against the state the host
 * gives them, and project windows and delta pages from it.
 */
export type FiberCommandKind = 'create' | 'append' | 'begin' | 'append-delta' | 'settle' | 'abort' | 'remove';

export type FiberKernelReceipt = Readonly<{
  readonly commandId: FiberCommandId;
  readonly fingerprint: string;
  readonly result: FiberCommandResult<unknown>;
}>;

export type FiberKernelFiber = FiberSnapshot &
  Readonly<{
    readonly limits: FiberLimits;
    readonly activeStreamCount: number;
  }>;

export type FiberKernelStream = Readonly<{
  readonly entryId: FiberEntryId;
  readonly nextOffset: number;
  readonly lastActivityAt: number;
  readonly retainedBytes: number;
  readonly retainedBatches: number;
  readonly status: 'active' | 'aborted-preview';
}>;

export type FiberKernelDecision<T> =
  | Readonly<{ readonly kind: 'result'; readonly persist: boolean; readonly result: FiberCommandResult<T> }>
  | Readonly<{ readonly kind: 'failure'; readonly failure: FiberFailure }>;

export type FiberKernelTouchDecision =
  | Readonly<{ readonly kind: 'result'; readonly persist: boolean; readonly result: FiberDeltaResult }>
  | Readonly<{ readonly kind: 'failure'; readonly failure: FiberFailure }>;

const CODEC_MAX_BYTES = 128;

function accepted<T>(commandId: FiberCommandId, value: T): FiberCommandResult<T> {
  return { status: 'accepted', commandId, value };
}

function failure(status: FiberFailure['status'], reason: string): FiberFailure {
  return { status, reason };
}

function replayOrConflict<T>(
  receipt: FiberKernelReceipt | undefined,
  fingerprint: string,
): FiberCommandResult<T> | FiberFailure | undefined {
  if (receipt === undefined) return undefined;
  if (receipt.fingerprint !== fingerprint) return failure('conflict', 'command-reused');
  return receipt.result as FiberCommandResult<T>;
}

export function fiberStoredEnvelopeBytes(envelope: FiberPayloadEnvelope): number {
  return new TextEncoder().encode(canonicalData(envelope)).byteLength;
}

export function validFiberPayloadEnvelope(value: unknown, maxBytes: number): value is FiberPayloadEnvelope {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) return false;
  const envelope = value as Readonly<Record<string, unknown>>;
  if (Object.keys(envelope).length !== 3) return false;
  if (typeof envelope.codec !== 'string' || envelope.codec.trim().length === 0) return false;
  if (new TextEncoder().encode(envelope.codec).byteLength > CODEC_MAX_BYTES) return false;
  if (envelope.version !== FIBER_PAYLOAD_VERSION) return false;
  try {
    assertStrictJsonData(envelope.value);
  } catch {
    return false;
  }
  return fiberStoredEnvelopeBytes(envelope as FiberPayloadEnvelope) <= maxBytes;
}

export function fiberCommandFingerprint(kind: FiberCommandKind, input: Readonly<Record<string, unknown>>): string {
  const { commandId: _commandId, ...rest } = input;
  return canonicalData({ kind, ...rest });
}

function activeFiber(fiber: FiberKernelFiber | undefined): FiberKernelFiber | FiberFailure {
  if (fiber === undefined) return failure('not-found', 'fiber-missing');
  if (fiber.lifecycle !== 'active') return failure('closed', 'fiber-closed');
  return fiber;
}

function snapshot(fiber: FiberKernelFiber): FiberSnapshot {
  return {
    fiberId: fiber.fiberId,
    ownerId: fiber.ownerId,
    nextPosition: fiber.nextPosition,
    lifecycle: fiber.lifecycle,
  };
}

function projectStreaming(entry: FiberEntry, stream: FiberKernelStream, now: number): FiberEntry {
  return {
    ...entry,
    state: { kind: 'streaming', nextOffset: stream.nextOffset, lastActivityAt: stream.lastActivityAt ?? now },
  };
}

export function prepareFiberCreate(input: {
  readonly command: CreateFiberCommand;
  readonly existingFiber?: FiberKernelFiber;
  readonly existingReceipt?: FiberKernelReceipt;
}): FiberKernelDecision<FiberSnapshot> {
  const fingerprint = fiberCommandFingerprint('create', {
    fiberId: input.command.fiberId,
    ownerId: input.command.ownerId,
    limits: input.command.limits,
  });
  const replay = replayOrConflict<FiberSnapshot>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as FiberFailure };
  }
  if (input.existingFiber !== undefined) return { kind: 'failure', failure: failure('conflict', 'fiber-exists') };
  let limits: FiberLimits;
  try {
    limits = createFiberLimits(input.command.limits, input.command.limits);
  } catch {
    return { kind: 'failure', failure: failure('invalid', 'limits') };
  }
  const fiber: FiberKernelFiber = {
    fiberId: input.command.fiberId,
    ownerId: input.command.ownerId,
    nextPosition: 0,
    lifecycle: 'active',
    limits,
    activeStreamCount: 0,
  };
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, snapshot(fiber)),
  };
}

export function prepareFiberAppendSettled(input: {
  readonly command: AppendSettledFiberEntryCommand;
  readonly now: number;
  readonly fiber?: FiberKernelFiber;
  readonly existingEntry?: FiberEntry;
  readonly existingReceipt?: FiberKernelReceipt;
}): FiberKernelDecision<FiberWriteResult> {
  const fingerprint = fiberCommandFingerprint('append', {
    fiberId: input.command.fiberId,
    entryId: input.command.entryId,
    descriptor: input.command.descriptor,
    payload: input.command.payload,
  });
  const replay = replayOrConflict<FiberWriteResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as FiberFailure };
  }
  const fiber = activeFiber(input.fiber);
  if ('status' in fiber) return { kind: 'failure', failure: fiber };
  if (!validFiberPayloadEnvelope(input.command.descriptor, fiber.limits.maxEntryBytes)) {
    return { kind: 'failure', failure: failure('invalid', 'descriptor') };
  }
  if (!validFiberPayloadEnvelope(input.command.payload, fiber.limits.maxEntryBytes)) {
    return { kind: 'failure', failure: failure('invalid', 'payload') };
  }
  if (input.existingEntry !== undefined) return { kind: 'failure', failure: failure('conflict', 'entry-exists') };
  const entry: FiberEntry = {
    entryId: input.command.entryId,
    fiberId: fiber.fiberId,
    position: fiber.nextPosition,
    descriptor: input.command.descriptor,
    state: { kind: 'settled', payload: input.command.payload, settledAt: input.now },
  };
  const next: FiberKernelFiber = { ...fiber, nextPosition: fiber.nextPosition + 1 };
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, { fiber: snapshot(next), entry }),
  };
}

export function prepareFiberBegin(input: {
  readonly command: BeginStreamingFiberEntryCommand;
  readonly now: number;
  readonly fiber?: FiberKernelFiber;
  readonly existingEntry?: FiberEntry;
  readonly existingReceipt?: FiberKernelReceipt;
}): FiberKernelDecision<FiberWriteResult> {
  const fingerprint = fiberCommandFingerprint('begin', {
    fiberId: input.command.fiberId,
    entryId: input.command.entryId,
    descriptor: input.command.descriptor,
  });
  const replay = replayOrConflict<FiberWriteResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as FiberFailure };
  }
  const fiber = activeFiber(input.fiber);
  if ('status' in fiber) return { kind: 'failure', failure: fiber };
  if (fiber.activeStreamCount >= fiber.limits.maxActiveStreams) {
    return { kind: 'failure', failure: failure('invalid', 'active-stream-limit') };
  }
  if (!validFiberPayloadEnvelope(input.command.descriptor, fiber.limits.maxEntryBytes)) {
    return { kind: 'failure', failure: failure('invalid', 'descriptor') };
  }
  if (input.existingEntry !== undefined) return { kind: 'failure', failure: failure('conflict', 'entry-exists') };
  const stream: FiberKernelStream = {
    entryId: input.command.entryId,
    nextOffset: 0,
    lastActivityAt: input.now,
    retainedBytes: 0,
    retainedBatches: 0,
    status: 'active',
  };
  const entry: FiberEntry = {
    entryId: input.command.entryId,
    fiberId: fiber.fiberId,
    position: fiber.nextPosition,
    descriptor: input.command.descriptor,
    state: { kind: 'streaming', nextOffset: 0, lastActivityAt: input.now },
  };
  const next: FiberKernelFiber = {
    ...fiber,
    nextPosition: fiber.nextPosition + 1,
    activeStreamCount: fiber.activeStreamCount + 1,
  };
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, {
      fiber: snapshot(next),
      entry: projectStreaming(entry, stream, input.now),
    }),
  };
}

export function prepareFiberAppendDelta(input: {
  readonly command: AppendFiberDeltaCommand;
  readonly now: number;
  readonly fiber?: FiberKernelFiber;
  readonly entry?: FiberEntry;
  readonly stream?: FiberKernelStream;
  readonly existingReceipt?: FiberKernelReceipt;
}): FiberKernelDecision<FiberAppendDeltaResult> {
  const fingerprint = fiberCommandFingerprint('append-delta', {
    fiberId: input.command.fiberId,
    entryId: input.command.entryId,
    start: input.command.start,
    payload: input.command.payload,
  });
  const replay = replayOrConflict<FiberAppendDeltaResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as FiberFailure };
  }
  const fiber = activeFiber(input.fiber);
  if ('status' in fiber) return { kind: 'failure', failure: fiber };
  if (input.entry === undefined || input.stream === undefined || input.stream.status !== 'active') {
    return { kind: 'failure', failure: failure(input.entry === undefined ? 'not-found' : 'invalid', 'not-streaming') };
  }
  if (input.entry.state.kind !== 'streaming') {
    return { kind: 'failure', failure: failure('invalid', 'terminal-entry') };
  }
  if (!validFiberPayloadEnvelope(input.command.payload, fiber.limits.maxBatchBytes)) {
    return { kind: 'failure', failure: failure('invalid', 'payload') };
  }
  if (input.command.start !== input.stream.nextOffset) {
    return { kind: 'failure', failure: failure('invalid', 'delta-offset') };
  }
  const batchBytes = fiberStoredEnvelopeBytes(input.command.payload);
  if (
    input.stream.retainedBytes + batchBytes > fiber.limits.maxActiveStreamBytes ||
    input.stream.retainedBatches + 1 > fiber.limits.maxActiveStreamBatches
  ) {
    const entry: FiberEntry = {
      ...input.entry,
      state: { kind: 'aborted', reason: 'limit-exceeded', endedAt: input.now },
    };
    const next: FiberKernelFiber = { ...fiber, activeStreamCount: Math.max(0, fiber.activeStreamCount - 1) };
    return {
      kind: 'result',
      persist: true,
      result: accepted(input.command.commandId, {
        kind: 'aborted',
        fiber: snapshot(next),
        entry,
        reason: 'limit-exceeded',
      }),
    };
  }
  const endExclusive = input.command.start + batchBytes;
  const stream: FiberKernelStream = {
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
      fiber: snapshot(fiber),
      entry,
      nextOffset: stream.nextOffset,
    }),
  };
}

export function prepareFiberHeartbeat(input: {
  readonly command: HeartbeatFiberEntryCommand;
  readonly now: number;
  readonly fiber?: FiberKernelFiber;
  readonly entry?: FiberEntry;
  readonly stream?: FiberKernelStream;
}): FiberKernelTouchDecision {
  const fiber = activeFiber(input.fiber);
  if ('status' in fiber) return { kind: 'failure', failure: fiber };
  if (input.entry === undefined || input.stream === undefined || input.stream.status !== 'active') {
    return { kind: 'failure', failure: failure(input.entry === undefined ? 'not-found' : 'invalid', 'not-streaming') };
  }
  const lastActivityAt = Math.max(input.stream.lastActivityAt, input.now);
  const stream: FiberKernelStream = { ...input.stream, lastActivityAt };
  const entry = projectStreaming(input.entry, stream, input.now);
  return {
    kind: 'result',
    persist: lastActivityAt > input.stream.lastActivityAt,
    result: {
      fiber: snapshot(fiber),
      entry,
      nextOffset: stream.nextOffset,
    },
  };
}

export function prepareFiberSettle(input: {
  readonly command: SettleFiberEntryCommand;
  readonly now: number;
  readonly fiber?: FiberKernelFiber;
  readonly entry?: FiberEntry;
  readonly stream?: FiberKernelStream;
  readonly existingReceipt?: FiberKernelReceipt;
}): FiberKernelDecision<FiberWriteResult> {
  const fingerprint = fiberCommandFingerprint('settle', {
    fiberId: input.command.fiberId,
    entryId: input.command.entryId,
    payload: input.command.payload,
  });
  const replay = replayOrConflict<FiberWriteResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as FiberFailure };
  }
  const fiber = activeFiber(input.fiber);
  if ('status' in fiber) return { kind: 'failure', failure: fiber };
  if (input.entry === undefined) return { kind: 'failure', failure: failure('not-found', 'entry-missing') };
  if (input.entry.state.kind === 'settled') {
    return { kind: 'failure', failure: failure('invalid', 'already-settled') };
  }
  if (input.entry.state.kind === 'aborted' || input.stream === undefined || input.stream.status !== 'active') {
    return { kind: 'failure', failure: failure('invalid', 'terminal-entry') };
  }
  if (!validFiberPayloadEnvelope(input.command.payload, fiber.limits.maxEntryBytes)) {
    return { kind: 'failure', failure: failure('invalid', 'payload') };
  }
  const entry: FiberEntry = {
    ...input.entry,
    state: { kind: 'settled', payload: input.command.payload, settledAt: input.now },
  };
  const next: FiberKernelFiber = { ...fiber, activeStreamCount: fiber.activeStreamCount - 1 };
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, { fiber: snapshot(next), entry }),
  };
}

export function prepareFiberAbort(input: {
  readonly command: AbortFiberEntryCommand;
  readonly now: number;
  readonly fiber?: FiberKernelFiber;
  readonly entry?: FiberEntry;
  readonly stream?: FiberKernelStream;
  readonly existingReceipt?: FiberKernelReceipt;
}): FiberKernelDecision<FiberWriteResult> {
  const fingerprint = fiberCommandFingerprint('abort', {
    fiberId: input.command.fiberId,
    entryId: input.command.entryId,
    reason: input.command.reason,
  });
  const replay = replayOrConflict<FiberWriteResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as FiberFailure };
  }
  const fiber = activeFiber(input.fiber);
  if ('status' in fiber) return { kind: 'failure', failure: fiber };
  if (input.entry === undefined) return { kind: 'failure', failure: failure('not-found', 'entry-missing') };
  if (input.entry.state.kind !== 'streaming' || input.stream === undefined || input.stream.status !== 'active') {
    return { kind: 'failure', failure: failure('invalid', 'terminal-entry') };
  }
  const entry: FiberEntry = {
    ...input.entry,
    state: { kind: 'aborted', reason: input.command.reason, endedAt: input.now },
  };
  const next: FiberKernelFiber = { ...fiber, activeStreamCount: fiber.activeStreamCount - 1 };
  return {
    kind: 'result',
    persist: true,
    result: accepted(input.command.commandId, { fiber: snapshot(next), entry }),
  };
}

export function prepareFiberRemove(input: {
  readonly command: RemoveFiberCommand;
  readonly fiber?: FiberKernelFiber;
  readonly existingReceipt?: FiberKernelReceipt;
}): FiberKernelDecision<FiberRemoveResult> {
  const fingerprint = fiberCommandFingerprint('remove', { fiberId: input.command.fiberId });
  const replay = replayOrConflict<FiberRemoveResult>(input.existingReceipt, fingerprint);
  if (replay !== undefined) {
    return 'status' in replay && replay.status === 'accepted'
      ? { kind: 'result', persist: false, result: replay }
      : { kind: 'failure', failure: replay as FiberFailure };
  }
  if (input.fiber === undefined) {
    return {
      kind: 'result',
      persist: false,
      result: accepted(input.command.commandId, { fiberId: input.command.fiberId, lifecycle: 'removed' }),
    };
  }
  switch (input.fiber.lifecycle) {
    case 'removed':
      return {
        kind: 'result',
        persist: false,
        result: accepted(input.command.commandId, { fiberId: input.command.fiberId, lifecycle: 'removed' }),
      };
    case 'closing':
      return {
        kind: 'result',
        persist: false,
        result: accepted(input.command.commandId, { fiberId: input.command.fiberId, lifecycle: 'closing' }),
      };
    case 'active':
      return {
        kind: 'result',
        persist: true,
        result: accepted(input.command.commandId, { fiberId: input.command.fiberId, lifecycle: 'closing' }),
      };
    default: {
      const exhaustive: never = input.fiber.lifecycle;
      throw new Error(`unexpected fiber lifecycle: ${exhaustive}`);
    }
  }
}

export function prepareFiberInactivityAbort(input: {
  readonly now: number;
  readonly fiber: FiberKernelFiber;
  readonly entry: FiberEntry;
  readonly stream: FiberKernelStream;
}):
  | Readonly<{ readonly kind: 'skip' }>
  | Readonly<{ readonly kind: 'abort'; readonly fiber: FiberSnapshot; readonly entry: FiberEntry }> {
  if (input.stream.status !== 'active' || input.entry.state.kind !== 'streaming') return { kind: 'skip' };
  if (input.now - input.stream.lastActivityAt < input.fiber.limits.inactivityTimeoutMs) return { kind: 'skip' };
  return {
    kind: 'abort',
    fiber: snapshot({ ...input.fiber, activeStreamCount: Math.max(0, input.fiber.activeStreamCount - 1) }),
    entry: {
      ...input.entry,
      state: { kind: 'aborted', reason: 'inactivity', endedAt: input.now },
    },
  };
}

function entryWindowBytes(entry: FiberEntry): number {
  const descriptor = fiberStoredEnvelopeBytes(entry.descriptor);
  return entry.state.kind === 'settled' ? descriptor + fiberStoredEnvelopeBytes(entry.state.payload) : descriptor;
}

export function fiberPageTakeBound(limit: number, ceiling: number): number {
  return Math.min(Math.max(limit, 0), ceiling) + 1;
}

export function fiberWindowRangeStart(
  request: Pick<ReadFiberWindowRequest, 'start' | 'direction'>,
  nextPosition: number,
): number {
  if (request.start !== undefined) return request.start;
  return request.direction === 'forward' ? 0 : nextPosition;
}

export function projectFiberWindow(input: {
  readonly fiber: FiberKernelFiber;
  readonly entries: readonly FiberEntry[];
  readonly streams: ReadonlyMap<string, FiberKernelStream>;
  readonly request: ReadFiberWindowRequest;
  readonly now: number;
}): FiberWindow | FiberFailure {
  if (input.request.limit < 1 || input.request.limit > input.fiber.limits.maxWindowEntries) {
    return failure('invalid', 'window-limit');
  }
  const ordered = [...input.entries].sort((left, right) => left.position - right.position);
  const start = fiberWindowRangeStart(input.request, input.fiber.nextPosition);
  const forward = input.request.direction === 'forward';
  const candidates = forward
    ? ordered.filter((entry) => entry.position >= start)
    : [...ordered.filter((entry) => entry.position < start)].reverse();
  const projected: FiberEntry[] = [];
  let bytes = 0;
  for (const entry of candidates) {
    const stream = input.streams.get(entry.entryId);
    const next =
      entry.state.kind === 'streaming' && stream !== undefined ? projectStreaming(entry, stream, input.now) : entry;
    const size = entryWindowBytes(next);
    if (projected.length >= input.request.limit) break;
    if (projected.length > 0 && bytes + size > input.fiber.limits.maxWindowBytes) break;
    if (projected.length === 0 && size > input.fiber.limits.maxWindowBytes) {
      return failure('invalid', 'window-bytes');
    }
    projected.push(next);
    bytes += size;
  }
  if (!forward) projected.reverse();
  const first = projected[0]?.position ?? start;
  const last = projected[projected.length - 1]?.position;
  return {
    fiberId: input.fiber.fiberId,
    entries: projected,
    start: first,
    endExclusive: last === undefined ? first : last + 1,
    hasMore: candidates.length > projected.length,
  };
}

export function projectFiberDeltas(input: {
  readonly fiber: FiberKernelFiber;
  readonly entry?: FiberEntry;
  readonly stream?: FiberKernelStream;
  readonly deltas: readonly FiberDeltaBatch[];
  readonly request: ReadFiberDeltasRequest;
}): FiberDeltaPage | FiberFailure | undefined {
  if (input.entry === undefined) return undefined;
  if (input.request.limit < 1 || input.request.limit > input.fiber.limits.maxDeltaPageBatches) {
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
  const page: FiberDeltaBatch[] = [];
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
    const size = fiberStoredEnvelopeBytes(batch.payload);
    if (page.length >= input.request.limit || bytes + size > input.fiber.limits.maxDeltaPageBytes) break;
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

export function nextFiberCleanupPage(input: {
  readonly cursor: FiberCleanupCursor;
  readonly ids: readonly string[];
  readonly pageSize: number;
}): FiberCleanupPage {
  const remaining = [...input.ids].filter((id) => input.cursor.after === undefined || id > input.cursor.after).sort();
  const ids = remaining.slice(0, input.pageSize);
  const last = ids[ids.length - 1];
  if (last !== undefined && remaining.length > ids.length) {
    return { kind: input.cursor.kind, ids, next: { kind: input.cursor.kind, after: last } };
  }
  return { kind: input.cursor.kind, ids };
}
