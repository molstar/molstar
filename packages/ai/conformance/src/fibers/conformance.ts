/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { FiberAppendDeltaResult, FiberClock, FiberCommandResult, FiberLimits } from '@molstar/ai-fibers';
import {
  type DeterministicFiberFixture,
  type FiberBackend,
  type FiberFixtureFactoryOptions,
  createDeterministicFiberFixture,
} from './fixtures.js';
import { fiberConformanceLimits } from './generators.js';

function accepted<T>(result: FiberCommandResult<T>): T {
  if (result.status !== 'accepted')
    throw new Error(`expected an accepted command, got ${result.status}: ${result.reason}`);
  return result.value;
}

function appended(
  result: FiberCommandResult<FiberAppendDeltaResult>,
): Extract<FiberAppendDeltaResult, { kind: 'appended' }> {
  const value = accepted(result);
  if (value.kind !== 'appended') throw new Error(`expected an appended delta, got ${value.kind}`);
  return value;
}

async function createdFiber(fixture: DeterministicFiberFixture, fiberId = fixture.fiberId): Promise<void> {
  accepted(
    await fixture.fibers.create({
      commandId: fixture.commandId(`create:${fiberId}`),
      fiberId,
      ownerId: fixture.ownerId,
      limits: fixture.limits,
    }),
  );
}

export const FIBER_CONFORMANCE_CASE_NAMES = [
  'replays exact Fiber commands and conflicts when input changes',
  'allocates monotonic positions for settled and streaming entries',
  'pages windows within Fiber limits',
  'requires contiguous deltas',
  'races settle and abort to one terminal',
  'returns settled windows without reconstructing deltas',
  'fences writes after Fiber removal',
  'aborts a stream that exceeds retained batch limits',
  'aborts a stream that goes silent past the inactivity timeout and frees its slot',
  'keeps entries, streams and receipts across a restart',
] as const;

/** What a backend gives the suite (analysis/viewer-mode.md §6). */
export interface FiberConformanceFixture {
  readonly name: string;
  /** A fiber port over fresh, empty storage for one case, on the case's clock and limits. */
  create(options: { readonly clock: FiberClock; readonly limits: FiberLimits }): FiberBackend;
  /** The storage outlives its port: `create`'s `reopen` gives a port over what the first one wrote. Default false. */
  readonly durable?: boolean;
}

/** Every Fiber conformance case, against one backend. Call it from a Jest file. */
export function runFiberConformance(backend: FiberConformanceFixture): void {
  const createFixture = (options: FiberFixtureFactoryOptions = {}): DeterministicFiberFixture =>
    createDeterministicFiberFixture({ createFiber: (deps) => backend.create(deps) }, options);

  describe(`Fiber conformance: ${backend.name}`, () => {
    it(FIBER_CONFORMANCE_CASE_NAMES[0], async () => {
      const fixture = createFixture();
      const input = {
        commandId: fixture.commandId('create'),
        fiberId: fixture.fiberId,
        ownerId: fixture.ownerId,
        limits: fixture.limits,
      };
      const first = await fixture.fibers.create(input);
      expect(await fixture.fibers.create(input)).toEqual(first);
      const conflict = await fixture.fibers.create({ ...input, fiberId: fixture.otherFiberId });
      expect(conflict.status).toBe('conflict');
    });

    it(FIBER_CONFORMANCE_CASE_NAMES[1], async () => {
      const fixture = createFixture();
      await createdFiber(fixture);
      const settled = accepted(
        await fixture.fibers.append({
          commandId: fixture.commandId('append'),
          fiberId: fixture.fiberId,
          entryId: fixture.entryId('e1'),
          descriptor: fixture.envelope({ role: 'left' }),
          payload: fixture.envelope('hello'),
        }),
      );
      const streaming = accepted(
        await fixture.fibers.begin({
          commandId: fixture.commandId('begin'),
          fiberId: fixture.fiberId,
          entryId: fixture.entryId('e2'),
          descriptor: fixture.envelope({ role: 'right' }),
        }),
      );
      expect(settled.entry.position).toBe(0);
      expect(streaming.entry.position).toBe(1);
      const window = await fixture.fibers.window({
        fiberId: fixture.fiberId,
        start: 0,
        limit: 8,
        direction: 'forward',
      });
      expect(window?.entries.map((entry) => entry.position)).toEqual([0, 1]);
    });

    it(FIBER_CONFORMANCE_CASE_NAMES[2], async () => {
      const fixture = createFixture();
      await createdFiber(fixture);
      await fixture.fibers.append({
        commandId: fixture.commandId('w1'),
        fiberId: fixture.fiberId,
        entryId: fixture.entryId('e1'),
        descriptor: fixture.envelope({ n: 1 }),
        payload: fixture.envelope('a'),
      });
      await fixture.fibers.append({
        commandId: fixture.commandId('w2'),
        fiberId: fixture.fiberId,
        entryId: fixture.entryId('e2'),
        descriptor: fixture.envelope({ n: 2 }),
        payload: fixture.envelope('b'),
      });
      const page = await fixture.fibers.window({
        fiberId: fixture.fiberId,
        start: 0,
        limit: 1,
        direction: 'forward',
      });
      expect(page?.entries.length).toBe(1);
      expect(page?.entries[0]?.position).toBe(0);
      expect(page?.hasMore).toBe(true);
    });

    it(FIBER_CONFORMANCE_CASE_NAMES[3], async () => {
      const fixture = createFixture();
      await createdFiber(fixture);
      const begun = accepted(
        await fixture.fibers.begin({
          commandId: fixture.commandId('begin-delta'),
          fiberId: fixture.fiberId,
          entryId: fixture.entryId('stream'),
          descriptor: fixture.envelope({ stream: true }),
        }),
      );
      const first = appended(
        await fixture.fibers.appendDelta({
          commandId: fixture.commandId('d1'),
          fiberId: fixture.fiberId,
          entryId: begun.entry.entryId,
          start: 0,
          payload: fixture.envelope('ab'),
        }),
      );
      const overlap = await fixture.fibers.appendDelta({
        commandId: fixture.commandId('d2'),
        fiberId: fixture.fiberId,
        entryId: begun.entry.entryId,
        start: 0,
        payload: fixture.envelope('cd'),
      });
      expect(overlap.status).toBe('invalid');
      const next = appended(
        await fixture.fibers.appendDelta({
          commandId: fixture.commandId('d3'),
          fiberId: fixture.fiberId,
          entryId: begun.entry.entryId,
          start: first.nextOffset,
          payload: fixture.envelope('ef'),
        }),
      );
      const page = await fixture.fibers.deltas({
        fiberId: fixture.fiberId,
        entryId: begun.entry.entryId,
        afterOffset: 0,
        limit: 1,
      });
      expect(page?.batches.length).toBe(1);
      expect(page?.hasMore).toBe(true);
      expect(next.entry.state.kind).toBe('streaming');
    });

    it(FIBER_CONFORMANCE_CASE_NAMES[4], async () => {
      const fixture = createFixture();
      await createdFiber(fixture);
      const begun = accepted(
        await fixture.fibers.begin({
          commandId: fixture.commandId('begin-race'),
          fiberId: fixture.fiberId,
          entryId: fixture.entryId('race'),
          descriptor: fixture.envelope({ stream: true }),
        }),
      );
      const settled = accepted(
        await fixture.fibers.settle({
          commandId: fixture.commandId('settle'),
          fiberId: fixture.fiberId,
          entryId: begun.entry.entryId,
          payload: fixture.envelope({ text: 'done' }),
        }),
      );
      expect(settled.entry.state.kind).toBe('settled');
      const abort = await fixture.fibers.abort({
        commandId: fixture.commandId('abort'),
        fiberId: fixture.fiberId,
        entryId: begun.entry.entryId,
        reason: 'caller-interrupted',
      });
      expect(abort.status).toBe('invalid');
      const late = await fixture.fibers.appendDelta({
        commandId: fixture.commandId('late-delta'),
        fiberId: fixture.fiberId,
        entryId: begun.entry.entryId,
        start: 0,
        payload: fixture.envelope('no'),
      });
      expect(late.status).not.toBe('accepted');
    });

    it(FIBER_CONFORMANCE_CASE_NAMES[5], async () => {
      const fixture = createFixture();
      await createdFiber(fixture);
      const begun = accepted(
        await fixture.fibers.begin({
          commandId: fixture.commandId('begin-settle-window'),
          fiberId: fixture.fiberId,
          entryId: fixture.entryId('preview'),
          descriptor: fixture.envelope({ stream: true }),
        }),
      );
      await fixture.fibers.appendDelta({
        commandId: fixture.commandId('preview-delta'),
        fiberId: fixture.fiberId,
        entryId: begun.entry.entryId,
        start: 0,
        payload: fixture.envelope('partial'),
      });
      const settled = accepted(
        await fixture.fibers.settle({
          commandId: fixture.commandId('settle-complete'),
          fiberId: fixture.fiberId,
          entryId: begun.entry.entryId,
          payload: fixture.envelope({ text: 'complete-from-caller' }),
        }),
      );
      const window = await fixture.fibers.window({
        fiberId: fixture.fiberId,
        start: 0,
        limit: 8,
        direction: 'forward',
      });
      expect(window?.entries[0]?.state.kind).toBe('settled');
      if (window?.entries[0]?.state.kind === 'settled') {
        expect(window.entries[0].state.payload).toEqual(fixture.envelope({ text: 'complete-from-caller' }));
      }
      const deltas = await fixture.fibers.deltas({
        fiberId: fixture.fiberId,
        entryId: begun.entry.entryId,
        afterOffset: 0,
        limit: 4,
      });
      expect(deltas?.batches.length ?? 0).toBe(0);
      expect(settled.entry.state.kind).toBe('settled');
    });

    it(FIBER_CONFORMANCE_CASE_NAMES[6], async () => {
      const fixture = createFixture();
      await createdFiber(fixture);
      const removed = accepted(
        await fixture.fibers.remove({ commandId: fixture.commandId('remove'), fiberId: fixture.fiberId }),
      );
      expect(removed.lifecycle === 'closing' || removed.lifecycle === 'removed').toBe(true);
      const late = await fixture.fibers.append({
        commandId: fixture.commandId('after-remove'),
        fiberId: fixture.fiberId,
        entryId: fixture.entryId('late'),
        descriptor: fixture.envelope({}),
        payload: fixture.envelope('no'),
      });
      expect(late.status).toBe('closed');
    });

    it(FIBER_CONFORMANCE_CASE_NAMES[7], async () => {
      const fixture = createFixture({ limits: fiberConformanceLimits({ maxActiveStreamBatches: 1 }) });
      await createdFiber(fixture);
      const begun = accepted(
        await fixture.fibers.begin({
          commandId: fixture.commandId('begin-limit'),
          fiberId: fixture.fiberId,
          entryId: fixture.entryId('limited'),
          descriptor: fixture.envelope({ stream: true }),
        }),
      );
      const first = appended(
        await fixture.fibers.appendDelta({
          commandId: fixture.commandId('limit-d1'),
          fiberId: fixture.fiberId,
          entryId: begun.entry.entryId,
          start: 0,
          payload: fixture.envelope('ab'),
        }),
      );
      const over = accepted(
        await fixture.fibers.appendDelta({
          commandId: fixture.commandId('limit-d2'),
          fiberId: fixture.fiberId,
          entryId: begun.entry.entryId,
          start: first.nextOffset,
          payload: fixture.envelope('cd'),
        }),
      );
      expect(over.kind).toBe('aborted');
      expect(over.kind === 'aborted' ? over.reason : undefined).toBe('limit-exceeded');
    });

    it(FIBER_CONFORMANCE_CASE_NAMES[8], async () => {
      const fixture = createFixture();
      await createdFiber(fixture);
      const { inactivityTimeoutMs, maxActiveStreams } = fixture.limits;
      const begin = async (name: string) =>
        accepted(
          await fixture.fibers.begin({
            commandId: fixture.commandId(`begin:${name}`),
            fiberId: fixture.fiberId,
            entryId: fixture.entryId(name),
            descriptor: fixture.envelope({ stream: name }),
          }),
        );
      const silent = await begin('silent');
      const beating = await begin('beating');
      appended(
        await fixture.fibers.appendDelta({
          commandId: fixture.commandId('silent-d1'),
          fiberId: fixture.fiberId,
          entryId: silent.entry.entryId,
          start: 0,
          payload: fixture.envelope('ab'),
        }),
      );
      // A heartbeat just before the timeout keeps its stream alive; the silent one is idle for the whole timeout.
      fixture.clock.advance(inactivityTimeoutMs - 1);
      const beat = await fixture.fibers.heartbeat({ fiberId: fixture.fiberId, entryId: beating.entry.entryId });
      expect('status' in beat).toBe(false);
      fixture.clock.advance(1);
      const window = await fixture.fibers.window({
        fiberId: fixture.fiberId,
        start: 0,
        limit: 8,
        direction: 'forward',
      });
      expect(window?.entries.map((entry) => entry.state.kind)).toEqual(['aborted', 'streaming']);
      const aborted = window?.entries[0]?.state;
      expect(aborted?.kind === 'aborted' ? aborted.reason : undefined).toBe('inactivity');
      // It is terminal: no late delta or settle is accepted.
      const late = await fixture.fibers.appendDelta({
        commandId: fixture.commandId('silent-late'),
        fiberId: fixture.fiberId,
        entryId: silent.entry.entryId,
        start: 0,
        payload: fixture.envelope('cd'),
      });
      expect(late.status).not.toBe('accepted');
      const settle = await fixture.fibers.settle({
        commandId: fixture.commandId('silent-settle'),
        fiberId: fixture.fiberId,
        entryId: silent.entry.entryId,
        payload: fixture.envelope({ text: 'too late' }),
      });
      expect(settle.status).toBe('invalid');
      // Its slot is free again: the live stream and `maxActiveStreams - 1` more fit, and one more does not.
      for (let n = 1; n < maxActiveStreams; n += 1) await begin(`extra-${n}`);
      const over = await fixture.fibers.begin({
        commandId: fixture.commandId('begin:over'),
        fiberId: fixture.fiberId,
        entryId: fixture.entryId('over'),
        descriptor: fixture.envelope({ stream: 'over' }),
      });
      expect(over).toMatchObject({ status: 'invalid', reason: 'active-stream-limit' });
    });

    (backend.durable === true ? it : it.skip)(FIBER_CONFORMANCE_CASE_NAMES[9], async () => {
      const fixture = createFixture();
      expect(fixture.reopen).toBeDefined();
      const create = {
        commandId: fixture.commandId('create'),
        fiberId: fixture.fiberId,
        ownerId: fixture.ownerId,
        limits: fixture.limits,
      };
      const created = accepted(await fixture.fibers.create(create));
      accepted(
        await fixture.fibers.append({
          commandId: fixture.commandId('append'),
          fiberId: fixture.fiberId,
          entryId: fixture.entryId('e1'),
          descriptor: fixture.envelope({ role: 'left' }),
          payload: fixture.envelope('hello'),
        }),
      );
      const begun = accepted(
        await fixture.fibers.begin({
          commandId: fixture.commandId('begin'),
          fiberId: fixture.fiberId,
          entryId: fixture.entryId('e2'),
          descriptor: fixture.envelope({ role: 'right' }),
        }),
      );
      const first = appended(
        await fixture.fibers.appendDelta({
          commandId: fixture.commandId('d1'),
          fiberId: fixture.fiberId,
          entryId: begun.entry.entryId,
          start: 0,
          payload: fixture.envelope('ab'),
        }),
      );
      const windowRequest = { fiberId: fixture.fiberId, start: 0, limit: 8, direction: 'forward' } as const;
      const deltaRequest = {
        fiberId: fixture.fiberId,
        entryId: begun.entry.entryId,
        afterOffset: 0,
        limit: 4,
      } as const;
      const window = await fixture.fibers.window(windowRequest);
      const deltas = await fixture.fibers.deltas(deltaRequest);
      expect(window?.entries.map((entry) => entry.state.kind)).toEqual(['settled', 'streaming']);
      expect(deltas?.batches.length).toBe(1);

      // A port over the same storage sees all of it, and answers a replayed command from the stored receipt.
      const restarted = fixture.reopen?.();
      if (restarted === undefined) throw new Error('a durable backend gives a reopen');
      expect(await restarted.window(windowRequest)).toEqual(window);
      expect(await restarted.deltas(deltaRequest)).toEqual(deltas);
      expect(await restarted.create(create)).toEqual({
        status: 'accepted',
        commandId: create.commandId,
        value: created,
      });
      // The stream goes on where it was, and settles.
      appended(
        await restarted.appendDelta({
          commandId: fixture.commandId('d2'),
          fiberId: fixture.fiberId,
          entryId: begun.entry.entryId,
          start: first.nextOffset,
          payload: fixture.envelope('cd'),
        }),
      );
      accepted(
        await restarted.settle({
          commandId: fixture.commandId('settle'),
          fiberId: fixture.fiberId,
          entryId: begun.entry.entryId,
          payload: fixture.envelope({ text: 'abcd' }),
        }),
      );
      // And so does that, for a third.
      const last = fixture.reopen?.();
      const settled = await last?.window(windowRequest);
      expect(settled?.entries.map((entry) => entry.state.kind)).toEqual(['settled', 'settled']);
      expect((await last?.deltas(deltaRequest))?.batches.length).toBe(0);
    });
  });
}
