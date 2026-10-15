/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';
import { type HarnessClient, HarnessError, type StopInput } from '../contracts/ports.js';
import {
  type HarnessEvent,
  type HarnessId,
  type HeadId,
  type LaunchKey,
  MAIN,
  type Seq,
  type StepEvent,
} from '../contracts/records.js';
import type { StorageTx } from '../contracts/storage.js';
import type { Host } from './index.js';
import type { Writer } from './writer.js';

type ClientOperations = Omit<HarnessClient, keyof import('../contracts/ports.js').HarnessReader>;

/** `{ seq, launch? }` from the events an operation wrote: the last seq, and its launch if any. */
export function launchResult(events: readonly HarnessEvent[]): { seq: Seq; launch?: LaunchKey } {
  const launched = events.find((event) => event.kind === 'control' && event.control === 'launched');
  const launch = launched?.kind === 'control' ? (launched.detail as { launch?: string }).launch : undefined;
  return { seq: events.at(-1)?.seq ?? 0, ...(launch !== undefined ? { launch } : {}) };
}

const done = (writer: Writer) => ({ seq: writer.exec.seq });
const withLaunch = (writer: Writer) => {
  const launch = writer.launchIfNeeded();
  return { seq: writer.exec.seq, ...(launch !== undefined ? { launch } : {}) };
};

/** Ends the active turn now, or asks its live runner to (harness.md §5). */
async function stopTurn(writer: Writer, reason: Json | undefined): Promise<void> {
  const turn = writer.turn();
  if (writer.liveLease()) {
    if (turn.stopRequested === undefined) {
      writer.set('turn', { ...turn, stopRequested: { at: writer.now, ...(reason !== undefined ? { reason } : {}) } });
      writer.control('stop-requested', { turnId: turn.id });
    }
    return;
  }
  await writer.endByControl('stopped', reason);
}

/** The head `main` had at the end of the last completed turn (reset's default target). */
async function lastCompletedMain(tx: StorageTx): Promise<HeadId | undefined> {
  let before: Seq | undefined;
  for (;;) {
    const page = await tx.events.endedTurns(before, 20);
    for (const event of page) {
      if (event.summary?.status !== 'completed' || event.kind !== 'step') continue;
      const steps = (await tx.events.byTurn(event.turnId, undefined, 10_000)).filter(
        (e): e is StepEvent => e.kind === 'step',
      );
      const moves = steps.flatMap((step) => step.moves.filter((move) => move.pointer === MAIN));
      return moves.at(-1)?.to;
    }
    if (page.length < 20) return undefined;
    before = page.at(-1)?.seq;
  }
}

export function clientOperations(host: Host, harnessId: HarnessId): ClientOperations {
  return {
    submit: (input) =>
      host.mutate(harnessId, {
        key: input.key,
        input: input as unknown as Json,
        replay: async (events) => launchResult(events),
        async run(writer) {
          writer.set('inbox', [
            ...writer.exec.inbox,
            {
              key: input.key,
              messages: input.messages,
              submittedAt: writer.now,
              ...(input.options !== undefined ? { options: input.options } : {}),
              ...(input.interrupt !== undefined ? { interrupt: input.interrupt } : {}),
            },
          ]);
          writer.control('submitted', { detail: { item: input.key } });
          return withLaunch(writer);
        },
      }),

    withdraw: (input) =>
      host.mutate(harnessId, {
        async run(writer) {
          if (!writer.exec.inbox.some((item) => item.key === input.key)) return { ...done(writer), withdrawn: false };
          writer.set(
            'inbox',
            writer.exec.inbox.filter((item) => item.key !== input.key),
          );
          writer.control('withdrawn', { detail: { item: input.key } });
          return { ...done(writer), withdrawn: true };
        },
      }),

    stop: (input: StopInput) =>
      host.mutate(harnessId, {
        key: input.key,
        input: input as unknown as Json,
        replay: async (events) => ({
          ...launchResult(events),
          ...(events.some((e) => e.kind === 'control' && e.control === 'withdrawn') ? { withdrawn: true } : {}),
        }),
        async run(writer) {
          if ('item' in input) {
            if (writer.exec.inbox.some((item) => item.key === input.item)) {
              writer.set(
                'inbox',
                writer.exec.inbox.filter((item) => item.key !== input.item),
              );
              writer.control('withdrawn', { detail: { item: input.item } });
              return { ...done(writer), withdrawn: true };
            }
            const turn = writer.exec.turn;
            if (turn === undefined) return done(writer);
            const events = await writer.tx.events.byTurn(turn.id, undefined, 10_000);
            const joined = events.some((e) => e.kind === 'step' && e.consumed.includes(input.item));
            if (!joined) return done(writer);
          } else if (writer.exec.turn?.id !== input.turnId) {
            throw new HarnessError('invalid', `turn ${input.turnId} is not active`);
          }
          await stopTurn(writer, input.reason);
          return withLaunch(writer);
        },
      }),

    resume: (input) =>
      host.mutate(harnessId, {
        key: input.key,
        input: input as unknown as Json,
        replay: async (events) => launchResult(events),
        async run(writer) {
          const turn = writer.exec.turn;
          if (turn?.id !== input.turnId || turn.park === undefined)
            throw new HarnessError('invalid', `turn ${input.turnId} is not parked`);
          const { park: _park, ...rest } = turn;
          writer.set('turn', rest);
          writer.control('resumed', { turnId: turn.id });
          return withLaunch(writer);
        },
      }),

    reset: (input) =>
      host.mutate(harnessId, {
        key: input.key,
        input: input as unknown as Json,
        replay: async (events) => launchResult(events),
        async run(writer) {
          const targets: Record<string, HeadId> = { ...input.pointers };
          if (input.pointers === undefined) {
            const main = await lastCompletedMain(writer.tx);
            if (main !== undefined) targets[MAIN] = main;
          }
          for (const [scope, head] of Object.entries(targets)) {
            const [record] = await writer.tx.heads.get([head]);
            if (record === undefined) throw new HarnessError('invalid', `no head ${head}`);
            if (scope === MAIN) writer.set('main', head);
            else writer.tx.pointers.put(scope, { scope, head });
          }
          if (writer.exec.turn !== undefined) await writer.endByControl('reset', input.reason);
          else
            writer.control('reset', {
              detail: { pointers: targets, ...(input.reason !== undefined ? { reason: input.reason } : {}) },
            });
          return withLaunch(writer);
        },
      }),

    check: () =>
      host.mutate(harnessId, {
        async run(writer) {
          const runner = writer.exec.runner;
          if (runner.phase === 'live' && runner.leaseUntil <= writer.now) {
            writer.set('runner', { phase: 'none' });
            writer.control('released', { detail: { runnerId: runner.runnerId, expired: true } });
          } else if (runner.phase === 'launching' && runner.until <= writer.now) {
            writer.set('runner', { phase: 'none' });
          }
          if (writer.exec.turn?.stopRequested !== undefined && !writer.liveLease())
            await writer.endByControl('stopped', writer.exec.turn.stopRequested.reason);
          return withLaunch(writer);
        },
      }),
  };
}
