/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { threadCommandId } from '../ids.js';
import {
  createThreadLimits,
  THREAD_CONFORMANCE_LIMITS,
  THREAD_LIMIT_HARD_BOUNDS,
  type ThreadAppendDeltaResult,
  type ThreadEntry,
  type ThreadKernelThread,
  type ThreadKernelReceipt,
  type ThreadKernelStream,
  type ThreadPayloadEnvelope,
  type ThreadSnapshot,
  threadCommandFingerprint,
  threadEntryId,
  threadId,
  threadOwnerId,
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
  validThreadPayloadEnvelope,
} from '../index.js';

const limits = createThreadLimits(THREAD_CONFORMANCE_LIMITS);
const owner = threadOwnerId('owner:thread');
const id = threadId('thread:1');
const envelope = (value: unknown): ThreadPayloadEnvelope => ({ codec: 'test.v1', version: 1, value });

function accepted<T>(decision: { kind: 'result' | 'failure'; result?: { status: string; value?: T } }): T {
  expect(decision.kind).toBe('result');
  if (decision.kind !== 'result' || decision.result?.status !== 'accepted') {
    throw new Error(`expected accepted, got ${JSON.stringify(decision)}`);
  }
  return decision.result.value as T;
}

function appended(decision: {
  kind: 'result' | 'failure';
  result?: { status: string; value?: ThreadAppendDeltaResult };
}): Extract<ThreadAppendDeltaResult, { kind: 'appended' }> {
  const value = accepted(decision);
  expect(value.kind).toBe('appended');
  if (value.kind !== 'appended') throw new Error('expected appended delta');
  return value;
}

function kernelThread(snapshot: ThreadSnapshot, activeStreamCount = 0): ThreadKernelThread {
  return { ...snapshot, limits, activeStreamCount };
}

describe('ThreadLimits', () => {
  it('accepts the conformance profile and rejects heartbeat after inactivity', () => {
    expect(createThreadLimits(THREAD_CONFORMANCE_LIMITS).heartbeatIntervalMs).toBe(50);
    expect(() =>
      createThreadLimits({ ...THREAD_CONFORMANCE_LIMITS, heartbeatIntervalMs: 200, inactivityTimeoutMs: 200 }),
    ).toThrow(/heartbeatIntervalMs/);
  });

  it('rejects a window that cannot hold one settled entry and accepts runtime ceilings', () => {
    expect(() => createThreadLimits({ ...THREAD_CONFORMANCE_LIMITS, maxWindowBytes: 4_096 })).toThrow(/maxWindowBytes/);
    expect(() =>
      createThreadLimits(THREAD_CONFORMANCE_LIMITS, { ...THREAD_LIMIT_HARD_BOUNDS, maxWindowEntries: 4 }),
    ).toThrow(/maxWindowEntries/);
    expect(
      createThreadLimits(
        { ...THREAD_CONFORMANCE_LIMITS, maxWindowEntries: 4 },
        { ...THREAD_LIMIT_HARD_BOUNDS, maxWindowEntries: 4 },
      ).maxWindowEntries,
    ).toBe(4);
  });
});

describe('Thread envelopes', () => {
  it('accepts version-1 strict JSON and rejects other shapes', () => {
    expect(validThreadPayloadEnvelope(envelope({ text: 'ok' }), 4_096)).toBe(true);
    expect(validThreadPayloadEnvelope({ codec: 'x', version: 2, value: null }, 4_096)).toBe(false);
    expect(validThreadPayloadEnvelope({ codec: 'x', version: 1, value: undefined }, 4_096)).toBe(false);
    expect(validThreadPayloadEnvelope({ codec: '', version: 1, value: null }, 4_096)).toBe(false);
  });
});

describe('Thread command fingerprint', () => {
  it('ignores command identity and conflicts when the rest changes', () => {
    const left = threadCommandFingerprint('append', { commandId: 'a', threadId: 'f', payload: envelope('x') });
    const right = threadCommandFingerprint('append', { commandId: 'b', threadId: 'f', payload: envelope('x') });
    const changed = threadCommandFingerprint('append', { commandId: 'a', threadId: 'f', payload: envelope('y') });
    expect(left).toBe(right);
    expect(left).not.toBe(changed);
  });
});

describe('Thread kernel', () => {
  it('replays exact create and conflicts when input changes', () => {
    const created = prepareThreadCreate({
      command: { commandId: threadCommandId('c1'), threadId: id, ownerId: owner, limits },
    });
    const snapshot = accepted<ThreadSnapshot>(created);
    expect(snapshot).toMatchObject({ threadId: id, nextPosition: 0, lifecycle: 'active' });
    const receipt: ThreadKernelReceipt = {
      commandId: threadCommandId('c1'),
      fingerprint: threadCommandFingerprint('create', { threadId: id, ownerId: owner, limits }),
      result: created.kind === 'result' ? created.result : { status: 'invalid', reason: 'x' },
    };
    const replay = prepareThreadCreate({
      command: { commandId: threadCommandId('c1'), threadId: id, ownerId: owner, limits },
      existingThread: kernelThread(snapshot),
      existingReceipt: receipt,
    });
    expect(replay).toMatchObject({ kind: 'result', persist: false });
    const conflict = prepareThreadCreate({
      command: { commandId: threadCommandId('c1'), threadId: threadId('thread:other'), ownerId: owner, limits },
      existingReceipt: receipt,
    });
    expect(conflict).toMatchObject({ kind: 'failure', failure: { status: 'conflict', reason: 'command-reused' } });
  });

  it('rejects invalid create limits before a Thread exists', () => {
    const invalid = prepareThreadCreate({
      command: {
        commandId: threadCommandId('bad-limits'),
        threadId: id,
        ownerId: owner,
        limits: { ...THREAD_CONFORMANCE_LIMITS, maxBatchBytes: THREAD_CONFORMANCE_LIMITS.maxEntryBytes + 1 },
      },
    });
    expect(invalid).toMatchObject({ kind: 'failure', failure: { status: 'invalid', reason: 'limits' } });
  });

  it('allocates monotonic positions for settled and streaming entries', () => {
    const thread = kernelThread(
      accepted(
        prepareThreadCreate({ command: { commandId: threadCommandId('c1'), threadId: id, ownerId: owner, limits } }),
      ),
    );
    const first = accepted(
      prepareThreadAppendSettled({
        command: {
          commandId: threadCommandId('c2'),
          threadId: id,
          entryId: threadEntryId('e1'),
          descriptor: envelope({ kind: 'user' }),
          payload: envelope({ text: 'hi' }),
        },
        now: 10,
        thread,
      }),
    );
    expect(first.entry.position).toBe(0);
    const streaming = accepted(
      prepareThreadBegin({
        command: {
          commandId: threadCommandId('c3'),
          threadId: id,
          entryId: threadEntryId('e2'),
          descriptor: envelope({ kind: 'assistant' }),
        },
        now: 11,
        thread: { ...thread, nextPosition: 1 },
      }),
    );
    expect(streaming.entry.position).toBe(1);
    expect(streaming.thread.nextPosition).toBe(2);
  });

  it('requires contiguous deltas and reports a gap on a hole', () => {
    const thread: ThreadKernelThread = {
      threadId: id,
      ownerId: owner,
      nextPosition: 1,
      lifecycle: 'active',
      limits,
      activeStreamCount: 1,
    };
    const entry: ThreadEntry = {
      entryId: threadEntryId('e2'),
      threadId: id,
      position: 0,
      descriptor: envelope({ kind: 'assistant' }),
      state: { kind: 'streaming', nextOffset: 0, lastActivityAt: 1 },
    };
    const stream: ThreadKernelStream = {
      entryId: entry.entryId,
      nextOffset: 0,
      lastActivityAt: 1,
      retainedBytes: 0,
      retainedBatches: 0,
      status: 'active',
    };
    const first = appended(
      prepareThreadAppendDelta({
        command: {
          commandId: threadCommandId('d1'),
          threadId: id,
          entryId: entry.entryId,
          start: 0,
          payload: envelope('ab'),
        },
        now: 2,
        thread,
        entry,
        stream,
      }),
    );
    expect(first.nextOffset).toBe(threadStoredEnvelopeBytes(envelope('ab')));
    const overlap = prepareThreadAppendDelta({
      command: {
        commandId: threadCommandId('d2'),
        threadId: id,
        entryId: entry.entryId,
        start: 0,
        payload: envelope('cd'),
      },
      now: 3,
      thread,
      entry: first.entry,
      stream: { ...stream, nextOffset: first.nextOffset, retainedBytes: 1, retainedBatches: 1 },
    });
    expect(overlap).toMatchObject({ kind: 'failure', failure: { reason: 'delta-offset' } });
    const page = projectThreadDeltas({
      thread,
      entry: first.entry,
      stream: { ...stream, nextOffset: first.nextOffset, status: 'active' },
      deltas: [
        { entryId: entry.entryId, start: 0, endExclusive: first.nextOffset, payload: envelope('ab') },
        {
          entryId: entry.entryId,
          start: first.nextOffset + 8,
          endExclusive: first.nextOffset + 16,
          payload: envelope('xx'),
        },
      ],
      request: { threadId: id, entryId: entry.entryId, afterOffset: first.nextOffset, limit: 4 },
    });
    expect(page).toMatchObject({ gap: { afterOffset: first.nextOffset } });
  });

  it('races settle and abort to one terminal and rejects late deltas', () => {
    const thread: ThreadKernelThread = {
      threadId: id,
      ownerId: owner,
      nextPosition: 1,
      lifecycle: 'active',
      limits,
      activeStreamCount: 1,
    };
    const entry: ThreadEntry = {
      entryId: threadEntryId('e2'),
      threadId: id,
      position: 0,
      descriptor: envelope({ kind: 'assistant' }),
      state: { kind: 'streaming', nextOffset: 4, lastActivityAt: 1 },
    };
    const stream: ThreadKernelStream = {
      entryId: entry.entryId,
      nextOffset: 4,
      lastActivityAt: 1,
      retainedBytes: 4,
      retainedBatches: 1,
      status: 'active',
    };
    const settled = accepted(
      prepareThreadSettle({
        command: {
          commandId: threadCommandId('s1'),
          threadId: id,
          entryId: entry.entryId,
          payload: envelope({ text: 'done' }),
        },
        now: 20,
        thread,
        entry,
        stream,
      }),
    );
    expect(settled.entry.state.kind).toBe('settled');
    const abort = prepareThreadAbort({
      command: { commandId: threadCommandId('a1'), threadId: id, entryId: entry.entryId, reason: 'caller-interrupted' },
      now: 21,
      thread: { ...thread, activeStreamCount: 0 },
      entry: settled.entry,
      stream,
    });
    expect(abort).toMatchObject({ kind: 'failure', failure: { reason: 'terminal-entry' } });
    const late = prepareThreadAppendDelta({
      command: {
        commandId: threadCommandId('d9'),
        threadId: id,
        entryId: entry.entryId,
        start: 4,
        payload: envelope('no'),
      },
      now: 22,
      thread: { ...thread, activeStreamCount: 0 },
      entry: settled.entry,
      stream,
    });
    expect(late).toMatchObject({ kind: 'failure' });
  });

  it('pages windows and fences remove', () => {
    const thread: ThreadKernelThread = {
      threadId: id,
      ownerId: owner,
      nextPosition: 2,
      lifecycle: 'active',
      limits,
      activeStreamCount: 0,
    };
    const entries: ThreadEntry[] = [
      {
        entryId: threadEntryId('e1'),
        threadId: id,
        position: 0,
        descriptor: envelope({ kind: 'user' }),
        state: { kind: 'settled', payload: envelope('a'), settledAt: 1 },
      },
      {
        entryId: threadEntryId('e2'),
        threadId: id,
        position: 1,
        descriptor: envelope({ kind: 'assistant' }),
        state: { kind: 'settled', payload: envelope('b'), settledAt: 2 },
      },
    ];
    const window = projectThreadWindow({
      thread,
      entries,
      streams: new Map(),
      request: { threadId: id, start: 0, limit: 1, direction: 'forward' },
      now: 3,
    });
    expect(window).toMatchObject({ hasMore: true, entries: [{ position: 0 }] });
    expect(threadPageTakeBound(1, limits.maxWindowEntries)).toBe(2);
    expect(threadWindowRangeStart({ direction: 'backward' }, 2)).toBe(2);
    const newest = projectThreadWindow({
      thread,
      entries,
      streams: new Map(),
      request: { threadId: id, limit: 1, direction: 'backward' },
      now: 3,
    });
    expect(newest).toMatchObject({ hasMore: true, entries: [{ position: 1 }] });
    const bounded = projectThreadWindow({
      thread: { ...thread, nextPosition: 256 },
      entries,
      streams: new Map(),
      request: { threadId: id, start: 0, limit: 1, direction: 'forward' },
      now: 3,
    });
    expect(bounded).toMatchObject({ hasMore: true, entries: [{ position: 0 }] });
    const closed = prepareThreadRemove({ command: { commandId: threadCommandId('r1'), threadId: id }, thread });
    expect(accepted(closed).lifecycle).toBe('closing');
    const fenced = prepareThreadAppendSettled({
      command: {
        commandId: threadCommandId('late'),
        threadId: id,
        entryId: threadEntryId('e9'),
        descriptor: envelope({}),
        payload: envelope('no'),
      },
      now: 4,
      thread: { ...thread, lifecycle: 'closing' },
    });
    expect(fenced).toMatchObject({ kind: 'failure', failure: { status: 'closed' } });
  });

  it('aborts a silent stream on inactivity and not sooner', () => {
    const thread: ThreadKernelThread = {
      threadId: id,
      ownerId: owner,
      nextPosition: 1,
      lifecycle: 'active',
      limits,
      activeStreamCount: 1,
    };
    const entry: ThreadEntry = {
      entryId: threadEntryId('e2'),
      threadId: id,
      position: 0,
      descriptor: envelope({ kind: 'assistant' }),
      state: { kind: 'streaming', nextOffset: 0, lastActivityAt: 0 },
    };
    const stream: ThreadKernelStream = {
      entryId: entry.entryId,
      nextOffset: 0,
      lastActivityAt: 0,
      retainedBytes: 0,
      retainedBatches: 0,
      status: 'active',
    };
    expect(prepareThreadInactivityAbort({ now: 199, thread, entry, stream })).toEqual({ kind: 'skip' });
    const aborted = prepareThreadInactivityAbort({ now: 200, thread, entry, stream });
    expect(aborted).toMatchObject({ kind: 'abort', entry: { state: { kind: 'aborted', reason: 'inactivity' } } });
    const beat = prepareThreadHeartbeat({
      command: { threadId: id, entryId: entry.entryId },
      now: 10,
      thread,
      entry,
      stream,
    });
    expect(beat).toMatchObject({
      kind: 'result',
      persist: true,
      result: { entry: { state: { kind: 'streaming', lastActivityAt: 10 } } },
    });
    const stale = prepareThreadHeartbeat({
      command: { threadId: id, entryId: entry.entryId },
      now: 4,
      thread,
      entry,
      stream: { ...stream, lastActivityAt: 10 },
    });
    expect(stale).toMatchObject({
      kind: 'result',
      persist: false,
      result: { entry: { state: { lastActivityAt: 10 } } },
    });
  });

  it('aborts a stream that would exceed retained batch or byte limits', () => {
    const streamLimits = createThreadLimits({ ...THREAD_CONFORMANCE_LIMITS, maxActiveStreamBatches: 1 });
    const thread: ThreadKernelThread = {
      threadId: id,
      ownerId: owner,
      nextPosition: 1,
      lifecycle: 'active',
      limits: streamLimits,
      activeStreamCount: 1,
    };
    const entry: ThreadEntry = {
      entryId: threadEntryId('e2'),
      threadId: id,
      position: 0,
      descriptor: envelope({ kind: 'assistant' }),
      state: { kind: 'streaming', nextOffset: 0, lastActivityAt: 1 },
    };
    const stream: ThreadKernelStream = {
      entryId: entry.entryId,
      nextOffset: 0,
      lastActivityAt: 1,
      retainedBytes: 0,
      retainedBatches: 1,
      status: 'active',
    };
    const aborted = accepted(
      prepareThreadAppendDelta({
        command: {
          commandId: threadCommandId('over'),
          threadId: id,
          entryId: entry.entryId,
          start: 0,
          payload: envelope('ab'),
        },
        now: 2,
        thread,
        entry,
        stream,
      }),
    );
    expect(aborted).toMatchObject({
      kind: 'aborted',
      reason: 'limit-exceeded',
      entry: { state: { kind: 'aborted', reason: 'limit-exceeded' } },
    });
    expect(aborted.kind === 'aborted' ? aborted.thread : undefined).toMatchObject({ lifecycle: 'active' });
  });
});
