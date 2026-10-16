/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  type BeginInput,
  type HarnessEvent,
  type LaunchKey,
  type Lease,
  MAIN,
  SCHEMA_VERSION,
} from '@molstar/ai-harness';
import { type ConformanceFixture, createSession, type Session } from '../session.js';
import { settle } from '../testing/clock.js';
import { user } from '../testing/messages.js';
import { appendMove, lastStep, TURN, with_ } from './helpers.js';

/** Polls (letting the runner and the simulated clock's promises settle) until the condition holds. */
async function until(condition: () => boolean, rounds = 200): Promise<void> {
  for (let i = 0; i < rounds; i += 1) {
    if (condition()) return;
    await settle();
  }
  throw new Error('until: the condition never held');
}

const statuses = async (s: Session): Promise<string[]> => (await s.turns()).map((turn) => turn.status);

/** A turn begun by hand through the runner port: its lease, ids and the head `main` has after `begin`. */
interface Started {
  readonly lease: Lease;
  readonly turnId: string;
  /** The pending step. */
  readonly step: string;
  readonly head: string;
  readonly begin: BeginInput;
}

type Wrap = <T>(call: () => Promise<T>) => Promise<T>;

/**
 * Submits `item-<tag>` (unless a launch is given), claims as `r-<tag>` and begins `hand-turn-<tag>`
 * with that item. `wrap` lets a case run each keyed call twice.
 */
async function startTurn(s: Session, tag: string, launch?: LaunchKey, wrap: Wrap = (call) => call()): Promise<Started> {
  const prompt = user(`hello ${tag}`);
  const given = launch ?? (await wrap(() => s.client.submit({ key: `item-${tag}`, messages: [prompt] }))).launch;
  if (given === undefined) throw new Error('startTurn: nothing launched');
  const { lease, execution } = await wrap(() =>
    s.port.claim({ launch: given, runnerId: `r-${tag}`, schemaVersion: SCHEMA_VERSION }),
  );
  const move = await appendMove(s.port, MAIN, execution.main, [prompt]);
  const begin: BeginInput = {
    key: `begin-${tag}`,
    turnId: `hand-turn-${tag}`,
    step: `hand-step-${tag}-1`,
    items: [`item-${tag}`],
    configuration: TURN,
    messages: [prompt],
    move,
    calls: [],
    next: { id: `hand-step-${tag}-2`, kind: 'model' },
  };
  const begun = await wrap(() => s.port.begin(lease, begin));
  return { lease: begun.lease, turnId: begin.turnId, step: `hand-step-${tag}-2`, head: move.head.id, begin };
}

/** T1–T6 (Turns and operations). */
export function turnCases(fixture: ConformanceFixture): void {
  describe('Turns and operations', () => {
    // T1. Operations run only under their preconditions (harness.md §12).

    it('T1 begin-while-active: begin with a turn active fails invalid', async () => {
      const s = await createSession(fixture);
      const started = await startTurn(s, '1');
      await s.client.submit({ key: 'item-2', messages: [user('hello 2')] });
      const move = await appendMove(s.port, MAIN, started.head, [user('hello 2')]);
      await expect(
        s.port.begin(started.lease, {
          ...started.begin,
          key: 'begin-2',
          turnId: 'hand-turn-2',
          step: 'hand-step-2-1',
          items: ['item-2'],
          messages: [user('hello 2')],
          move,
          next: { id: 'hand-step-2-2', kind: 'model' },
        }),
      ).rejects.toMatchObject({ code: 'invalid' });
    });

    // T4. Stop ends the turn and touches nothing else (harness.md §5).

    it("T4 stop-live-runner: Stop with a live runner commits abandoned with the ended attempts' usage and end stopped", async () => {
      const s = await createSession(fixture, { model: [{ text: 'slow', waitMs: 60_000 }] });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);

      await s.client.stop({ key: 'stop-1', turnId: (await s.client.execution()).turn?.id as string });
      await s.clock.runUntil(running);

      const last = lastStep(await s.events());
      expect(last).toMatchObject({ ending: 'abandoned', end: { status: 'stopped' } });
      const attempts = last.calls.flatMap((call) => call.attempts);
      expect(attempts).toHaveLength(1);
      expect(attempts[0]).toMatchObject({ outcome: 'aborted' });
      expect(attempts[0]?.usage).toBeDefined();
      expect(last.summary?.usage.modelCalls).toBeGreaterThanOrEqual(1);
      expect(await statuses(s)).toEqual(['stopped']);
    });

    it("T4 stop-beats-settlement: Stop requested before a settlement's commit refuses the commit and abandons the step", async () => {
      const s: Session = await createSession(fixture, {
        model: [{ text: 'final' }],
        drive: with_({
          // The last decision before the settling commit: Stop lands just before it.
          boundary: async () => {
            await s.client.stop({ key: 'stop-1', turnId: (await s.client.execution()).turn?.id as string });
            return 'proceed' as const;
          },
        }),
      });
      await s.send('go');

      expect(await statuses(s)).toEqual(['stopped']);
      expect(lastStep(await s.events())).toMatchObject({ ending: 'abandoned', end: { status: 'stopped' } });
      expect((await s.main()).map((message) => message.role)).toEqual(['user']);
    });

    // T5. Reset ends any turn and moves only the heads it names (harness.md §5).

    /** A completed turn, then a second one begun by hand, whose `begin` moved `main` past the first. */
    const completedThenBegun = async (s: Session) => {
      await s.send('first');
      const completed = (await s.client.execution()).main;
      const completedMessages = await s.main();
      await startTurn(s, '2');
      const during = (await s.client.execution()).main;
      const duringMessages = await s.main();
      expect(during).not.toBe(completed);
      return { completed, completedMessages, during, duringMessages };
    };

    it('T5 reset-default: reset without heads moves main to the head at the end of the last completed turn', async () => {
      const s = await createSession(fixture, { model: [{ text: 'one' }] });
      const { completed, completedMessages } = await completedThenBegun(s);

      await s.client.reset({ key: 'reset-1' });

      expect((await s.client.execution()).main).toBe(completed);
      expect(await s.main()).toEqual(completedMessages);
      expect(await statuses(s)).toEqual(['reset', 'completed']);
    });

    // T6. The event log is complete and ordered (harness.md §3, §12, §13).

    it('T6 catch-up: a subscriber joining mid-turn subscribes, then reads events(after), and sees each event once', async () => {
      let clock: { sleep(ms: number, signal?: AbortSignal): Promise<void> } | undefined;
      const s = await createSession(fixture, {
        tools: {
          slow: async (context) => {
            await clock?.sleep(5_000, context.signal);
            return { text: 'ok' };
          },
        },
        model: [{ calls: [{ tool: 'slow', arguments: {} }] }, { text: 'final' }],
      });
      clock = s.clock;
      const { launch } = await s.send('go', { run: false });
      const runner = s.runner();
      const running = runner.run(launch as string);
      await until(() => s.tools.callsOf('slow').length === 1);

      const seen: number[] = [];
      const observed: number[] = [];
      let last = 0;
      const take = (event: HarnessEvent) => {
        if (event.seq <= last) return;
        last = event.seq;
        seen.push(event.seq);
      };
      runner.observe().subscribe((observation) => {
        if (observation.kind !== 'event') return;
        observed.push(observation.event.seq);
        take(observation.event);
      });
      for (const event of await s.client.events({ after: last })) take(event);
      await s.clock.runUntil(running);
      for (const event of await s.client.events({ after: last })) take(event);

      expect(observed.length).toBeGreaterThan(0);
      expect(seen).toEqual((await s.events()).map((event) => event.seq));
    });
  });
}
