/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { HarnessEvent } from '@molstar/ai-harness';
import { createSession, type Session } from '../session.js';
import { lastStep, until, with_ } from './helpers.js';

const statuses = async (s: Session): Promise<string[]> => (await s.turns()).map((turn) => turn.status);

/** T1–T6 (Turns and operations). */
export function turnCases(): void {
  describe('Turns and operations', () => {
    // T1. Operations run only under their preconditions (harness.md §5).

    it('T1 begin-while-active: a turn while a turn is active is refused invalid', async () => {
      const s = createSession({ model: [{ text: 'slow', waitMs: 60_000 }] });
      const { handle } = await s.send('go', { wait: false });
      expect(() => s.harness.turn({ messages: [] })).toThrow(expect.objectContaining({ code: 'invalid' }));
      await until(() => s.model.invocations.length === 1);
      expect(() => s.harness.turn({ messages: [] })).toThrow(expect.objectContaining({ code: 'invalid' }));
      await s.run(handle);
      expect(await statuses(s)).toEqual(['completed']);
    });

    // T4. Stop ends the turn and touches nothing else (harness.md §5).

    it("T4 stop-live-runner: Stop during a step records it abandoned with the ended attempts' usage and end stopped", async () => {
      const s = createSession({ model: [{ text: 'slow', waitMs: 60_000 }] });
      const { handle } = await s.send('go', { wait: false });
      await until(() => s.model.invocations.length === 1);

      s.harness.stop({ by: 'person' });
      expect(await s.run(handle)).toEqual({
        turnId: handle.turnId,
        end: { status: 'stopped', reason: { by: 'person' } },
      });

      const last = lastStep(await s.events());
      expect(last).toMatchObject({ ending: 'abandoned', end: { status: 'stopped' } });
      const attempts = last.calls.flatMap((call) => call.attempts);
      expect(attempts).toHaveLength(1);
      expect(attempts[0]).toMatchObject({ outcome: 'aborted' });
      expect(attempts[0]?.usage).toBeDefined();
      expect(last.summary?.usage.modelCalls).toBeGreaterThanOrEqual(1);
      expect(await statuses(s)).toEqual(['stopped']);
    });

    it('T4 stop-beats-settlement: Stop asked before a settlement is written abandons the step', async () => {
      const s: Session = createSession({
        model: [{ text: 'final' }],
        drive: with_({
          // The last decision before the settling write: Stop lands just before it.
          boundary: () => {
            s.harness.stop();
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

    it('T5 reset-default: reset without heads moves main to the head at the end of the last completed turn', async () => {
      const s = createSession({ model: [{ text: 'one' }, { text: 'slow', waitMs: 60_000 }] });
      await s.send('first');
      const completed = (await s.harness.execution()).main;
      const completedMessages = await s.main();
      const { handle } = await s.send('second', { wait: false });
      await until(() => s.model.invocations.length === 2);
      expect((await s.harness.execution()).main).not.toBe(completed);

      s.harness.reset();

      expect(await s.run(handle)).toEqual({ turnId: handle.turnId, end: { status: 'reset' } });
      expect((await s.harness.execution()).main).toBe(completed);
      expect(await s.main()).toEqual(completedMessages);
      expect(await statuses(s)).toEqual(['reset', 'completed']);
    });

    // T6. The event log is complete and ordered (harness.md §3, §13).

    it('T6 catch-up: a subscriber joining mid-turn subscribes, then reads events(after), and sees each event once', async () => {
      let clock: { sleep(ms: number, signal?: AbortSignal): Promise<void> } | undefined;
      const s = createSession({
        tools: {
          slow: async (context) => {
            await clock?.sleep(5_000, context.signal);
            return { text: 'ok' };
          },
        },
        model: [{ calls: [{ tool: 'slow', arguments: {} }] }, { text: 'final' }],
      });
      clock = s.clock;
      const { handle } = await s.send('go', { wait: false });
      await until(() => s.tools.callsOf('slow').length === 1);

      const seen: number[] = [];
      const observed: number[] = [];
      let last = 0;
      const take = (event: HarnessEvent) => {
        if (event.seq <= last) return;
        last = event.seq;
        seen.push(event.seq);
      };
      s.harness.observe().subscribe((observation) => {
        if (observation.kind !== 'event') return;
        observed.push(observation.event.seq);
        take(observation.event);
      });
      for (const event of await s.harness.events({ after: last })) take(event);
      await s.run(handle);
      for (const event of await s.harness.events({ after: last })) take(event);

      expect(observed.length).toBeGreaterThan(0);
      expect(seen).toEqual((await s.events()).map((event) => event.seq));
    });
  });
}
