/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { fiberCommandId } from '../ids.js';
import {
  createFiberLimits,
  FIBER_CONFORMANCE_LIMITS,
  FIBER_LIMIT_HARD_BOUNDS,
  type FiberAppendDeltaResult,
  type FiberEntry,
  type FiberKernelFiber,
  type FiberKernelReceipt,
  type FiberKernelStream,
  type FiberPayloadEnvelope,
  type FiberSnapshot,
  fiberCommandFingerprint,
  fiberEntryId,
  fiberId,
  fiberOwnerId,
  fiberPageTakeBound,
  fiberStoredEnvelopeBytes,
  fiberWindowRangeStart,
  nextFiberCleanupPage,
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
  validFiberPayloadEnvelope,
} from '../index.js';

const limits = createFiberLimits(FIBER_CONFORMANCE_LIMITS);
const owner = fiberOwnerId('owner:fiber');
const id = fiberId('fiber:1');
const envelope = (value: unknown): FiberPayloadEnvelope => ({ codec: 'test.v1', version: 1, value });

function accepted<T>(decision: { kind: 'result' | 'failure'; result?: { status: string; value?: T } }): T {
  expect(decision.kind).toBe('result');
  if (decision.kind !== 'result' || decision.result?.status !== 'accepted') {
    throw new Error(`expected accepted, got ${JSON.stringify(decision)}`);
  }
  return decision.result.value as T;
}

function appended(decision: {
  kind: 'result' | 'failure';
  result?: { status: string; value?: FiberAppendDeltaResult };
}): Extract<FiberAppendDeltaResult, { kind: 'appended' }> {
  const value = accepted(decision);
  expect(value.kind).toBe('appended');
  if (value.kind !== 'appended') throw new Error('expected appended delta');
  return value;
}

function kernelFiber(snapshot: FiberSnapshot, activeStreamCount = 0): FiberKernelFiber {
  return { ...snapshot, limits, activeStreamCount };
}

describe('FiberLimits', () => {
  it('accepts the conformance profile and rejects heartbeat after inactivity', () => {
    expect(createFiberLimits(FIBER_CONFORMANCE_LIMITS).heartbeatIntervalMs).toBe(50);
    expect(() =>
      createFiberLimits({ ...FIBER_CONFORMANCE_LIMITS, heartbeatIntervalMs: 200, inactivityTimeoutMs: 200 }),
    ).toThrow(/heartbeatIntervalMs/);
  });

  it('rejects a window that cannot hold one settled entry and accepts runtime ceilings', () => {
    expect(() => createFiberLimits({ ...FIBER_CONFORMANCE_LIMITS, maxWindowBytes: 4_096 })).toThrow(/maxWindowBytes/);
    expect(() =>
      createFiberLimits(FIBER_CONFORMANCE_LIMITS, { ...FIBER_LIMIT_HARD_BOUNDS, maxWindowEntries: 4 }),
    ).toThrow(/maxWindowEntries/);
    expect(
      createFiberLimits(
        { ...FIBER_CONFORMANCE_LIMITS, maxWindowEntries: 4 },
        { ...FIBER_LIMIT_HARD_BOUNDS, maxWindowEntries: 4 },
      ).maxWindowEntries,
    ).toBe(4);
  });
});

describe('Fiber envelopes', () => {
  it('accepts version-1 strict JSON and rejects other shapes', () => {
    expect(validFiberPayloadEnvelope(envelope({ text: 'ok' }), 4_096)).toBe(true);
    expect(validFiberPayloadEnvelope({ codec: 'x', version: 2, value: null }, 4_096)).toBe(false);
    expect(validFiberPayloadEnvelope({ codec: 'x', version: 1, value: undefined }, 4_096)).toBe(false);
    expect(validFiberPayloadEnvelope({ codec: '', version: 1, value: null }, 4_096)).toBe(false);
  });
});

describe('Fiber command fingerprint', () => {
  it('ignores command identity and conflicts when the rest changes', () => {
    const left = fiberCommandFingerprint('append', { commandId: 'a', fiberId: 'f', payload: envelope('x') });
    const right = fiberCommandFingerprint('append', { commandId: 'b', fiberId: 'f', payload: envelope('x') });
    const changed = fiberCommandFingerprint('append', { commandId: 'a', fiberId: 'f', payload: envelope('y') });
    expect(left).toBe(right);
    expect(left).not.toBe(changed);
  });
});

describe('Fiber kernel', () => {
  it('replays exact create and conflicts when input changes', () => {
    const created = prepareFiberCreate({
      command: { commandId: fiberCommandId('c1'), fiberId: id, ownerId: owner, limits },
    });
    const snapshot = accepted<FiberSnapshot>(created);
    expect(snapshot).toMatchObject({ fiberId: id, nextPosition: 0, lifecycle: 'active' });
    const receipt: FiberKernelReceipt = {
      commandId: fiberCommandId('c1'),
      fingerprint: fiberCommandFingerprint('create', { fiberId: id, ownerId: owner, limits }),
      result: created.kind === 'result' ? created.result : { status: 'invalid', reason: 'x' },
    };
    const replay = prepareFiberCreate({
      command: { commandId: fiberCommandId('c1'), fiberId: id, ownerId: owner, limits },
      existingFiber: kernelFiber(snapshot),
      existingReceipt: receipt,
    });
    expect(replay).toMatchObject({ kind: 'result', persist: false });
    const conflict = prepareFiberCreate({
      command: { commandId: fiberCommandId('c1'), fiberId: fiberId('fiber:other'), ownerId: owner, limits },
      existingReceipt: receipt,
    });
    expect(conflict).toMatchObject({ kind: 'failure', failure: { status: 'conflict', reason: 'command-reused' } });
  });

  it('rejects invalid create limits before a Fiber exists', () => {
    const invalid = prepareFiberCreate({
      command: {
        commandId: fiberCommandId('bad-limits'),
        fiberId: id,
        ownerId: owner,
        limits: { ...FIBER_CONFORMANCE_LIMITS, maxBatchBytes: FIBER_CONFORMANCE_LIMITS.maxEntryBytes + 1 },
      },
    });
    expect(invalid).toMatchObject({ kind: 'failure', failure: { status: 'invalid', reason: 'limits' } });
  });

  it('allocates monotonic positions for settled and streaming entries', () => {
    const fiber = kernelFiber(
      accepted(
        prepareFiberCreate({ command: { commandId: fiberCommandId('c1'), fiberId: id, ownerId: owner, limits } }),
      ),
    );
    const first = accepted(
      prepareFiberAppendSettled({
        command: {
          commandId: fiberCommandId('c2'),
          fiberId: id,
          entryId: fiberEntryId('e1'),
          descriptor: envelope({ kind: 'user' }),
          payload: envelope({ text: 'hi' }),
        },
        now: 10,
        fiber,
      }),
    );
    expect(first.entry.position).toBe(0);
    const streaming = accepted(
      prepareFiberBegin({
        command: {
          commandId: fiberCommandId('c3'),
          fiberId: id,
          entryId: fiberEntryId('e2'),
          descriptor: envelope({ kind: 'assistant' }),
        },
        now: 11,
        fiber: { ...fiber, nextPosition: 1 },
      }),
    );
    expect(streaming.entry.position).toBe(1);
    expect(streaming.fiber.nextPosition).toBe(2);
  });

  it('requires contiguous deltas and reports a gap on a hole', () => {
    const fiber: FiberKernelFiber = {
      fiberId: id,
      ownerId: owner,
      nextPosition: 1,
      lifecycle: 'active',
      limits,
      activeStreamCount: 1,
    };
    const entry: FiberEntry = {
      entryId: fiberEntryId('e2'),
      fiberId: id,
      position: 0,
      descriptor: envelope({ kind: 'assistant' }),
      state: { kind: 'streaming', nextOffset: 0, lastActivityAt: 1 },
    };
    const stream: FiberKernelStream = {
      entryId: entry.entryId,
      nextOffset: 0,
      lastActivityAt: 1,
      retainedBytes: 0,
      retainedBatches: 0,
      status: 'active',
    };
    const first = appended(
      prepareFiberAppendDelta({
        command: {
          commandId: fiberCommandId('d1'),
          fiberId: id,
          entryId: entry.entryId,
          start: 0,
          payload: envelope('ab'),
        },
        now: 2,
        fiber,
        entry,
        stream,
      }),
    );
    expect(first.nextOffset).toBe(fiberStoredEnvelopeBytes(envelope('ab')));
    const overlap = prepareFiberAppendDelta({
      command: {
        commandId: fiberCommandId('d2'),
        fiberId: id,
        entryId: entry.entryId,
        start: 0,
        payload: envelope('cd'),
      },
      now: 3,
      fiber,
      entry: first.entry,
      stream: { ...stream, nextOffset: first.nextOffset, retainedBytes: 1, retainedBatches: 1 },
    });
    expect(overlap).toMatchObject({ kind: 'failure', failure: { reason: 'delta-offset' } });
    const page = projectFiberDeltas({
      fiber,
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
      request: { fiberId: id, entryId: entry.entryId, afterOffset: first.nextOffset, limit: 4 },
    });
    expect(page).toMatchObject({ gap: { afterOffset: first.nextOffset } });
  });

  it('races settle and abort to one terminal and rejects late deltas', () => {
    const fiber: FiberKernelFiber = {
      fiberId: id,
      ownerId: owner,
      nextPosition: 1,
      lifecycle: 'active',
      limits,
      activeStreamCount: 1,
    };
    const entry: FiberEntry = {
      entryId: fiberEntryId('e2'),
      fiberId: id,
      position: 0,
      descriptor: envelope({ kind: 'assistant' }),
      state: { kind: 'streaming', nextOffset: 4, lastActivityAt: 1 },
    };
    const stream: FiberKernelStream = {
      entryId: entry.entryId,
      nextOffset: 4,
      lastActivityAt: 1,
      retainedBytes: 4,
      retainedBatches: 1,
      status: 'active',
    };
    const settled = accepted(
      prepareFiberSettle({
        command: {
          commandId: fiberCommandId('s1'),
          fiberId: id,
          entryId: entry.entryId,
          payload: envelope({ text: 'done' }),
        },
        now: 20,
        fiber,
        entry,
        stream,
      }),
    );
    expect(settled.entry.state.kind).toBe('settled');
    const abort = prepareFiberAbort({
      command: { commandId: fiberCommandId('a1'), fiberId: id, entryId: entry.entryId, reason: 'caller-interrupted' },
      now: 21,
      fiber: { ...fiber, activeStreamCount: 0 },
      entry: settled.entry,
      stream,
    });
    expect(abort).toMatchObject({ kind: 'failure', failure: { reason: 'terminal-entry' } });
    const late = prepareFiberAppendDelta({
      command: {
        commandId: fiberCommandId('d9'),
        fiberId: id,
        entryId: entry.entryId,
        start: 4,
        payload: envelope('no'),
      },
      now: 22,
      fiber: { ...fiber, activeStreamCount: 0 },
      entry: settled.entry,
      stream,
    });
    expect(late).toMatchObject({ kind: 'failure' });
  });

  it('pages windows, cleanup cursors, and fences remove', () => {
    const fiber: FiberKernelFiber = {
      fiberId: id,
      ownerId: owner,
      nextPosition: 2,
      lifecycle: 'active',
      limits,
      activeStreamCount: 0,
    };
    const entries: FiberEntry[] = [
      {
        entryId: fiberEntryId('e1'),
        fiberId: id,
        position: 0,
        descriptor: envelope({ kind: 'user' }),
        state: { kind: 'settled', payload: envelope('a'), settledAt: 1 },
      },
      {
        entryId: fiberEntryId('e2'),
        fiberId: id,
        position: 1,
        descriptor: envelope({ kind: 'assistant' }),
        state: { kind: 'settled', payload: envelope('b'), settledAt: 2 },
      },
    ];
    const window = projectFiberWindow({
      fiber,
      entries,
      streams: new Map(),
      request: { fiberId: id, start: 0, limit: 1, direction: 'forward' },
      now: 3,
    });
    expect(window).toMatchObject({ hasMore: true, entries: [{ position: 0 }] });
    expect(fiberPageTakeBound(1, limits.maxWindowEntries)).toBe(2);
    expect(fiberWindowRangeStart({ direction: 'backward' }, 2)).toBe(2);
    const newest = projectFiberWindow({
      fiber,
      entries,
      streams: new Map(),
      request: { fiberId: id, limit: 1, direction: 'backward' },
      now: 3,
    });
    expect(newest).toMatchObject({ hasMore: true, entries: [{ position: 1 }] });
    const bounded = projectFiberWindow({
      fiber: { ...fiber, nextPosition: 256 },
      entries,
      streams: new Map(),
      request: { fiberId: id, start: 0, limit: 1, direction: 'forward' },
      now: 3,
    });
    expect(bounded).toMatchObject({ hasMore: true, entries: [{ position: 0 }] });
    const page = nextFiberCleanupPage({
      cursor: { kind: 'entries' },
      ids: ['e1', 'e2', 'e3'],
      pageSize: 2,
    });
    expect(page.ids).toEqual(['e1', 'e2']);
    expect(page.next).toEqual({ kind: 'entries', after: 'e2' });
    const closed = prepareFiberRemove({ command: { commandId: fiberCommandId('r1'), fiberId: id }, fiber });
    expect(accepted(closed).lifecycle).toBe('closing');
    const fenced = prepareFiberAppendSettled({
      command: {
        commandId: fiberCommandId('late'),
        fiberId: id,
        entryId: fiberEntryId('e9'),
        descriptor: envelope({}),
        payload: envelope('no'),
      },
      now: 4,
      fiber: { ...fiber, lifecycle: 'closing' },
    });
    expect(fenced).toMatchObject({ kind: 'failure', failure: { status: 'closed' } });
  });

  it('aborts a silent stream on inactivity and not sooner', () => {
    const fiber: FiberKernelFiber = {
      fiberId: id,
      ownerId: owner,
      nextPosition: 1,
      lifecycle: 'active',
      limits,
      activeStreamCount: 1,
    };
    const entry: FiberEntry = {
      entryId: fiberEntryId('e2'),
      fiberId: id,
      position: 0,
      descriptor: envelope({ kind: 'assistant' }),
      state: { kind: 'streaming', nextOffset: 0, lastActivityAt: 0 },
    };
    const stream: FiberKernelStream = {
      entryId: entry.entryId,
      nextOffset: 0,
      lastActivityAt: 0,
      retainedBytes: 0,
      retainedBatches: 0,
      status: 'active',
    };
    expect(prepareFiberInactivityAbort({ now: 199, fiber, entry, stream })).toEqual({ kind: 'skip' });
    const aborted = prepareFiberInactivityAbort({ now: 200, fiber, entry, stream });
    expect(aborted).toMatchObject({ kind: 'abort', entry: { state: { kind: 'aborted', reason: 'inactivity' } } });
    const beat = prepareFiberHeartbeat({
      command: { fiberId: id, entryId: entry.entryId },
      now: 10,
      fiber,
      entry,
      stream,
    });
    expect(beat).toMatchObject({
      kind: 'result',
      persist: true,
      result: { entry: { state: { kind: 'streaming', lastActivityAt: 10 } } },
    });
    const stale = prepareFiberHeartbeat({
      command: { fiberId: id, entryId: entry.entryId },
      now: 4,
      fiber,
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
    const streamLimits = createFiberLimits({ ...FIBER_CONFORMANCE_LIMITS, maxActiveStreamBatches: 1 });
    const fiber: FiberKernelFiber = {
      fiberId: id,
      ownerId: owner,
      nextPosition: 1,
      lifecycle: 'active',
      limits: streamLimits,
      activeStreamCount: 1,
    };
    const entry: FiberEntry = {
      entryId: fiberEntryId('e2'),
      fiberId: id,
      position: 0,
      descriptor: envelope({ kind: 'assistant' }),
      state: { kind: 'streaming', nextOffset: 0, lastActivityAt: 1 },
    };
    const stream: FiberKernelStream = {
      entryId: entry.entryId,
      nextOffset: 0,
      lastActivityAt: 1,
      retainedBytes: 0,
      retainedBatches: 1,
      status: 'active',
    };
    const aborted = accepted(
      prepareFiberAppendDelta({
        command: {
          commandId: fiberCommandId('over'),
          fiberId: id,
          entryId: entry.entryId,
          start: 0,
          payload: envelope('ab'),
        },
        now: 2,
        fiber,
        entry,
        stream,
      }),
    );
    expect(aborted).toMatchObject({
      kind: 'aborted',
      reason: 'limit-exceeded',
      entry: { state: { kind: 'aborted', reason: 'limit-exceeded' } },
    });
    expect(aborted.kind === 'aborted' ? aborted.fiber : undefined).toMatchObject({ lifecycle: 'active' });
  });
});
