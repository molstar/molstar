/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Clock } from '../../messages/index.js';
import {
  THREAD_LIMIT_HARD_BOUNDS,
  type ThreadAppendDeltaResult,
  type ThreadCommandResult,
  type ThreadFailure,
  type ThreadLimits,
} from '../../threads/index.js';
import {
  type DeterministicThreadFixture,
  type ThreadBackend,
  type ThreadFixtureFactoryOptions,
  createDeterministicThreadFixture,
} from './fixtures.js';
import { threadConformanceLimits } from './generators.js';

/** A read's page; a failure fails the case, naming why. */
function read<T extends object>(result: T | ThreadFailure): T {
  if ('status' in result && 'reason' in result) throw new Error(`read failed: ${result.status} ${result.reason}`);
  return result as T;
}

function accepted<T>(result: ThreadCommandResult<T>): T {
  if (result.status !== 'accepted')
    throw new Error(`expected an accepted command, got ${result.status}: ${result.reason}`);
  return result.value;
}

function appended(
  result: ThreadCommandResult<ThreadAppendDeltaResult>,
): Extract<ThreadAppendDeltaResult, { kind: 'appended' }> {
  const value = accepted(result);
  if (value.kind !== 'appended') throw new Error(`expected an appended delta, got ${value.kind}`);
  return value;
}

async function createdThread(fixture: DeterministicThreadFixture, threadId = fixture.threadId): Promise<void> {
  accepted(
    await fixture.threads.create({
      commandId: fixture.commandId(`create:${threadId}`),
      threadId,
      ownerId: fixture.ownerId,
      limits: fixture.limits,
    }),
  );
}

export const THREAD_CONFORMANCE_CASE_NAMES = [
  'replays exact Thread commands and conflicts when input changes',
  'allocates monotonic positions for settled and streaming entries',
  'pages windows within Thread limits',
  'requires contiguous deltas',
  'races settle and abort to one terminal',
  'returns settled windows without reconstructing deltas',
  'fences writes after Thread removal',
  'aborts a stream that exceeds retained batch limits',
  'aborts a stream that goes silent past the inactivity timeout and frees its slot',
  'keeps entries, streams and receipts across a restart',
  'bounds the limits a Thread is created with',
  'says why a read fails: a missing Thread or entry, a request out of its limits',
] as const;

/** What a backend gives the suite (threads.md § Storage). */
export interface ThreadConformanceFixture {
  readonly name: string;
  /** A thread port over fresh, empty storage for one case, on the case's clock and limits. */
  create(options: { readonly clock: Pick<Clock, 'now'>; readonly limits: ThreadLimits }): ThreadBackend;
  /** The storage outlives its port: `create`'s `reopen` gives a port over what the first one wrote. Default false. */
  readonly durable?: boolean;
}

/** Every Thread conformance case, against one backend. Call it from a Vitest file. */
export function runThreadConformance(backend: ThreadConformanceFixture): void {
  const createFixture = (options: ThreadFixtureFactoryOptions = {}): DeterministicThreadFixture =>
    createDeterministicThreadFixture({ createThread: (deps) => backend.create(deps) }, options);

  describe(`Thread conformance: ${backend.name}`, () => {
    it(THREAD_CONFORMANCE_CASE_NAMES[0], async () => {
      const fixture = createFixture();
      const input = {
        commandId: fixture.commandId('create'),
        threadId: fixture.threadId,
        ownerId: fixture.ownerId,
        limits: fixture.limits,
      };
      const first = await fixture.threads.create(input);
      expect(await fixture.threads.create(input)).toEqual(first);
      const conflict = await fixture.threads.create({ ...input, threadId: fixture.otherThreadId });
      expect(conflict.status).toBe('conflict');
    });

    it(THREAD_CONFORMANCE_CASE_NAMES[1], async () => {
      const fixture = createFixture();
      await createdThread(fixture);
      const settled = accepted(
        await fixture.threads.append({
          commandId: fixture.commandId('append'),
          threadId: fixture.threadId,
          entryId: fixture.entryId('e1'),
          descriptor: fixture.envelope({ role: 'left' }),
          payload: fixture.envelope('hello'),
        }),
      );
      const streaming = accepted(
        await fixture.threads.begin({
          commandId: fixture.commandId('begin'),
          threadId: fixture.threadId,
          entryId: fixture.entryId('e2'),
          descriptor: fixture.envelope({ role: 'right' }),
        }),
      );
      expect(settled.entry.position).toBe(0);
      expect(streaming.entry.position).toBe(1);
      const window = read(
        await fixture.threads.window({
          threadId: fixture.threadId,
          start: 0,
          limit: 8,
          direction: 'forward',
        }),
      );
      expect(window?.entries.map((entry) => entry.position)).toEqual([0, 1]);
    });

    it(THREAD_CONFORMANCE_CASE_NAMES[2], async () => {
      const fixture = createFixture();
      await createdThread(fixture);
      await fixture.threads.append({
        commandId: fixture.commandId('w1'),
        threadId: fixture.threadId,
        entryId: fixture.entryId('e1'),
        descriptor: fixture.envelope({ n: 1 }),
        payload: fixture.envelope('a'),
      });
      await fixture.threads.append({
        commandId: fixture.commandId('w2'),
        threadId: fixture.threadId,
        entryId: fixture.entryId('e2'),
        descriptor: fixture.envelope({ n: 2 }),
        payload: fixture.envelope('b'),
      });
      const page = read(
        await fixture.threads.window({
          threadId: fixture.threadId,
          start: 0,
          limit: 1,
          direction: 'forward',
        }),
      );
      expect(page?.entries.length).toBe(1);
      expect(page?.entries[0]?.position).toBe(0);
      expect(page?.hasMore).toBe(true);
    });

    it(THREAD_CONFORMANCE_CASE_NAMES[3], async () => {
      const fixture = createFixture();
      await createdThread(fixture);
      const begun = accepted(
        await fixture.threads.begin({
          commandId: fixture.commandId('begin-delta'),
          threadId: fixture.threadId,
          entryId: fixture.entryId('stream'),
          descriptor: fixture.envelope({ stream: true }),
        }),
      );
      const first = appended(
        await fixture.threads.appendDelta({
          commandId: fixture.commandId('d1'),
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          start: 0,
          payload: fixture.envelope('ab'),
        }),
      );
      const overlap = await fixture.threads.appendDelta({
        commandId: fixture.commandId('d2'),
        threadId: fixture.threadId,
        entryId: begun.entry.entryId,
        start: 0,
        payload: fixture.envelope('cd'),
      });
      expect(overlap.status).toBe('invalid');
      const next = appended(
        await fixture.threads.appendDelta({
          commandId: fixture.commandId('d3'),
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          start: first.nextOffset,
          payload: fixture.envelope('ef'),
        }),
      );
      const page = read(
        await fixture.threads.deltas({
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          afterOffset: 0,
          limit: 1,
        }),
      );
      expect(page?.batches.length).toBe(1);
      expect(page?.hasMore).toBe(true);
      expect(next.entry.state.kind).toBe('streaming');
    });

    it(THREAD_CONFORMANCE_CASE_NAMES[4], async () => {
      const fixture = createFixture();
      await createdThread(fixture);
      const begun = accepted(
        await fixture.threads.begin({
          commandId: fixture.commandId('begin-race'),
          threadId: fixture.threadId,
          entryId: fixture.entryId('race'),
          descriptor: fixture.envelope({ stream: true }),
        }),
      );
      const settled = accepted(
        await fixture.threads.settle({
          commandId: fixture.commandId('settle'),
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          payload: fixture.envelope({ text: 'done' }),
        }),
      );
      expect(settled.entry.state.kind).toBe('settled');
      const abort = await fixture.threads.abort({
        commandId: fixture.commandId('abort'),
        threadId: fixture.threadId,
        entryId: begun.entry.entryId,
        reason: 'caller-interrupted',
      });
      expect(abort.status).toBe('invalid');
      const late = await fixture.threads.appendDelta({
        commandId: fixture.commandId('late-delta'),
        threadId: fixture.threadId,
        entryId: begun.entry.entryId,
        start: 0,
        payload: fixture.envelope('no'),
      });
      expect(late.status).not.toBe('accepted');
    });

    it(THREAD_CONFORMANCE_CASE_NAMES[5], async () => {
      const fixture = createFixture();
      await createdThread(fixture);
      const begun = accepted(
        await fixture.threads.begin({
          commandId: fixture.commandId('begin-settle-window'),
          threadId: fixture.threadId,
          entryId: fixture.entryId('preview'),
          descriptor: fixture.envelope({ stream: true }),
        }),
      );
      await fixture.threads.appendDelta({
        commandId: fixture.commandId('preview-delta'),
        threadId: fixture.threadId,
        entryId: begun.entry.entryId,
        start: 0,
        payload: fixture.envelope('partial'),
      });
      const settled = accepted(
        await fixture.threads.settle({
          commandId: fixture.commandId('settle-complete'),
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          payload: fixture.envelope({ text: 'complete-from-caller' }),
        }),
      );
      const window = read(
        await fixture.threads.window({
          threadId: fixture.threadId,
          start: 0,
          limit: 8,
          direction: 'forward',
        }),
      );
      expect(window?.entries[0]?.state.kind).toBe('settled');
      if (window?.entries[0]?.state.kind === 'settled') {
        expect(window.entries[0].state.payload).toEqual(fixture.envelope({ text: 'complete-from-caller' }));
      }
      const deltas = read(
        await fixture.threads.deltas({
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          afterOffset: 0,
          limit: 4,
        }),
      );
      expect(deltas?.batches.length ?? 0).toBe(0);
      expect(settled.entry.state.kind).toBe('settled');
    });

    it(THREAD_CONFORMANCE_CASE_NAMES[6], async () => {
      const fixture = createFixture();
      await createdThread(fixture);
      const removed = accepted(
        await fixture.threads.remove({ commandId: fixture.commandId('remove'), threadId: fixture.threadId }),
      );
      expect(removed.lifecycle === 'closing' || removed.lifecycle === 'removed').toBe(true);
      const late = await fixture.threads.append({
        commandId: fixture.commandId('after-remove'),
        threadId: fixture.threadId,
        entryId: fixture.entryId('late'),
        descriptor: fixture.envelope({}),
        payload: fixture.envelope('no'),
      });
      expect(late.status).toBe('closed');
    });

    it(THREAD_CONFORMANCE_CASE_NAMES[7], async () => {
      const fixture = createFixture({ limits: threadConformanceLimits({ maxActiveStreamBatches: 1 }) });
      await createdThread(fixture);
      const begun = accepted(
        await fixture.threads.begin({
          commandId: fixture.commandId('begin-limit'),
          threadId: fixture.threadId,
          entryId: fixture.entryId('limited'),
          descriptor: fixture.envelope({ stream: true }),
        }),
      );
      const first = appended(
        await fixture.threads.appendDelta({
          commandId: fixture.commandId('limit-d1'),
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          start: 0,
          payload: fixture.envelope('ab'),
        }),
      );
      const over = accepted(
        await fixture.threads.appendDelta({
          commandId: fixture.commandId('limit-d2'),
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          start: first.nextOffset,
          payload: fixture.envelope('cd'),
        }),
      );
      expect(over.kind).toBe('aborted');
      expect(over.kind === 'aborted' ? over.reason : undefined).toBe('limit-exceeded');
    });

    it(THREAD_CONFORMANCE_CASE_NAMES[8], async () => {
      const fixture = createFixture();
      await createdThread(fixture);
      const { inactivityTimeoutMs, maxActiveStreams } = fixture.limits;
      const begin = async (name: string) =>
        accepted(
          await fixture.threads.begin({
            commandId: fixture.commandId(`begin:${name}`),
            threadId: fixture.threadId,
            entryId: fixture.entryId(name),
            descriptor: fixture.envelope({ stream: name }),
          }),
        );
      const silent = await begin('silent');
      const beating = await begin('beating');
      appended(
        await fixture.threads.appendDelta({
          commandId: fixture.commandId('silent-d1'),
          threadId: fixture.threadId,
          entryId: silent.entry.entryId,
          start: 0,
          payload: fixture.envelope('ab'),
        }),
      );
      // A heartbeat just before the timeout keeps its stream alive; the silent one is idle for the whole timeout.
      fixture.clock.advance(inactivityTimeoutMs - 1);
      const beat = await fixture.threads.heartbeat({ threadId: fixture.threadId, entryId: beating.entry.entryId });
      expect('status' in beat).toBe(false);
      fixture.clock.advance(1);
      const window = read(
        await fixture.threads.window({
          threadId: fixture.threadId,
          start: 0,
          limit: 8,
          direction: 'forward',
        }),
      );
      expect(window?.entries.map((entry) => entry.state.kind)).toEqual(['aborted', 'streaming']);
      const aborted = window?.entries[0]?.state;
      expect(aborted?.kind === 'aborted' ? aborted.reason : undefined).toBe('inactivity');
      // It is terminal: no late delta or settle is accepted.
      const late = await fixture.threads.appendDelta({
        commandId: fixture.commandId('silent-late'),
        threadId: fixture.threadId,
        entryId: silent.entry.entryId,
        start: 0,
        payload: fixture.envelope('cd'),
      });
      expect(late.status).not.toBe('accepted');
      const settle = await fixture.threads.settle({
        commandId: fixture.commandId('silent-settle'),
        threadId: fixture.threadId,
        entryId: silent.entry.entryId,
        payload: fixture.envelope({ text: 'too late' }),
      });
      expect(settle.status).toBe('invalid');
      // Its slot is free again: the live stream and `maxActiveStreams - 1` more fit, and one more does not.
      for (let n = 1; n < maxActiveStreams; n += 1) await begin(`extra-${n}`);
      const over = await fixture.threads.begin({
        commandId: fixture.commandId('begin:over'),
        threadId: fixture.threadId,
        entryId: fixture.entryId('over'),
        descriptor: fixture.envelope({ stream: 'over' }),
      });
      expect(over).toMatchObject({ status: 'invalid', reason: 'active-stream-limit' });
    });

    (backend.durable === true ? it : it.skip)(THREAD_CONFORMANCE_CASE_NAMES[9], async () => {
      const fixture = createFixture();
      expect(fixture.reopen).toBeDefined();
      const create = {
        commandId: fixture.commandId('create'),
        threadId: fixture.threadId,
        ownerId: fixture.ownerId,
        limits: fixture.limits,
      };
      const created = accepted(await fixture.threads.create(create));
      accepted(
        await fixture.threads.append({
          commandId: fixture.commandId('append'),
          threadId: fixture.threadId,
          entryId: fixture.entryId('e1'),
          descriptor: fixture.envelope({ role: 'left' }),
          payload: fixture.envelope('hello'),
        }),
      );
      const begun = accepted(
        await fixture.threads.begin({
          commandId: fixture.commandId('begin'),
          threadId: fixture.threadId,
          entryId: fixture.entryId('e2'),
          descriptor: fixture.envelope({ role: 'right' }),
        }),
      );
      const first = appended(
        await fixture.threads.appendDelta({
          commandId: fixture.commandId('d1'),
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          start: 0,
          payload: fixture.envelope('ab'),
        }),
      );
      const windowRequest = { threadId: fixture.threadId, start: 0, limit: 8, direction: 'forward' } as const;
      const deltaRequest = {
        threadId: fixture.threadId,
        entryId: begun.entry.entryId,
        afterOffset: 0,
        limit: 4,
      } as const;
      const window = read(await fixture.threads.window(windowRequest));
      const deltas = read(await fixture.threads.deltas(deltaRequest));
      expect(window?.entries.map((entry) => entry.state.kind)).toEqual(['settled', 'streaming']);
      expect(deltas?.batches.length).toBe(1);

      // A port over the same storage sees all of it, and answers a replayed command from the stored receipt.
      const restarted = fixture.reopen?.();
      if (restarted === undefined) throw new Error('a durable backend gives a reopen');
      expect(read(await restarted.window(windowRequest))).toEqual(window);
      expect(read(await restarted.deltas(deltaRequest))).toEqual(deltas);
      expect(await restarted.create(create)).toEqual({
        status: 'accepted',
        commandId: create.commandId,
        value: created,
      });
      // The stream goes on where it was, and settles.
      appended(
        await restarted.appendDelta({
          commandId: fixture.commandId('d2'),
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          start: first.nextOffset,
          payload: fixture.envelope('cd'),
        }),
      );
      accepted(
        await restarted.settle({
          commandId: fixture.commandId('settle'),
          threadId: fixture.threadId,
          entryId: begun.entry.entryId,
          payload: fixture.envelope({ text: 'abcd' }),
        }),
      );
      // And so does that, for a third.
      const last = fixture.reopen?.();
      const settled = last && read(await last.window(windowRequest));
      expect(settled?.entries.map((entry) => entry.state.kind)).toEqual(['settled', 'settled']);
      expect((last && read(await last.deltas(deltaRequest)))?.batches.length).toBe(0);
    });

    it(THREAD_CONFORMANCE_CASE_NAMES[10], async () => {
      const fixture = createFixture();
      const over = { ...fixture.limits, maxActiveStreams: THREAD_LIMIT_HARD_BOUNDS.maxActiveStreams + 1 };
      const refused = await fixture.threads.create({
        commandId: fixture.commandId('create'),
        threadId: fixture.threadId,
        ownerId: fixture.ownerId,
        limits: over,
      });
      expect(refused).toEqual({ status: 'invalid', reason: 'limits' });
    });

    it(THREAD_CONFORMANCE_CASE_NAMES[11], async () => {
      const fixture = createFixture();
      const window = { threadId: fixture.threadId, limit: 1, direction: 'forward' } as const;
      expect(await fixture.threads.window(window)).toEqual({ status: 'not-found', reason: 'thread-missing' });
      await createdThread(fixture);
      const entryId = fixture.entryId('e1');
      expect(await fixture.threads.deltas({ threadId: fixture.threadId, entryId, afterOffset: 0, limit: 1 })).toEqual({
        status: 'not-found',
        reason: 'entry-missing',
      });
      expect(await fixture.threads.window({ ...window, limit: fixture.limits.maxWindowEntries + 1 })).toEqual({
        status: 'invalid',
        reason: 'window-limit',
      });
      accepted(
        await fixture.threads.begin({
          commandId: fixture.commandId('begin'),
          threadId: fixture.threadId,
          entryId,
          descriptor: fixture.envelope({ role: 'reply' }),
        }),
      );
      const deltas = { threadId: fixture.threadId, entryId, afterOffset: 0, limit: 0 };
      expect(await fixture.threads.deltas(deltas)).toEqual({ status: 'invalid', reason: 'delta-page-limit' });
    });
  });
}
