/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { LaunchKey, RunnerPort } from '@molstar/ai-harness';
import { type ConformanceFixture, createSession, type Session } from '../session.js';
import { settle } from '../testing/clock.js';
import { user } from '../testing/messages.js';
import { delegate, stepEvents, with_ } from './helpers.js';

/** Handlers run inside the session they are made for; this is their way back to it. */
interface Later {
  s?: Session;
}

async function until(condition: () => boolean, limit = 500): Promise<void> {
  for (let i = 0; i < limit && !condition(); i += 1) await settle();
  if (!condition()) throw new Error('until: the condition never held');
}

const runWith = (s: Session, launch: LaunchKey | undefined, port: RunnerPort) =>
  s.run(launch as string, s.runner({ port }));

const scripted = { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} };

/** R1–R12: the phase-4 review's proven findings (progress.md), each kept as a case. */
export function regressionCases(fixture: ConformanceFixture): void {
  describe('Regressions', () => {
    it('R1 cold-read-fault: a failed head-record read is retried, not a loss (harness.md §4)', async () => {
      const s = await createSession(fixture, { model: [{ text: 'ok' }] });
      let failed = false;
      const port = delegate(s.port, {
        headRecords: async (ids) => {
          if (!failed) {
            failed = true;
            throw new Error('read failed');
          }
          return await s.port.headRecords(ids);
        },
      });
      const { launch } = await s.send('go', { run: false });
      expect((await runWith(s, launch, port)).reason).toBe('released');
      expect((await s.turns())[0]?.status).toBe('completed');
    });

    it('R2 interrupt-already-stored: an interrupt whose text is already in main commits interrupted (harness.md §10)', async () => {
      const s = await createSession(fixture, {
        model: [{ text: 'slow', waitMs: 60_000 }, { text: 'after' }],
      });
      const { launch } = await s.send('wait', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);
      const turnId = (await s.client.execution()).turn?.id as string;
      await s.client.submit({ key: 'int-1', messages: [user('wait')], interrupt: { turnId } });
      expect((await s.clock.runUntil(running)).reason).toBe('released');
      expect(stepEvents(await s.events()).some((e) => e.ending === 'interrupted')).toBe(true);
      expect((await s.turns())[0]?.status).toBe('completed');
    });

    it("R3 drive-run-in-model-result: a drive's nested run does not replace the step's message (harness.md §6, §7)", async () => {
      const s = await createSession(fixture, {
        tools: { a: () => ({ text: 'a done' }) },
        model: [
          { calls: [{ tool: 'a', arguments: {} }], finish: 'output-limit' },
          { text: 'nested summary' },
          { text: 'final' },
        ],
        drive: with_({
          'model-result': async (_d, context) => {
            if (context.depth === 0)
              await context.run({
                input: [user('summarize')],
                tools: { enabled: [], preauthorized: [], definitions: [] },
                model: scripted,
                limits: {},
              });
            return 'accept' as const;
          },
        }),
      });
      await s.send('go');
      const withCall = (await s.main()).filter(
        (m) => m.role === 'assistant' && m.parts.some((p) => p.kind === 'tool-call'),
      );
      expect(withCall).toHaveLength(1);
    });

    it("R4 stop-refuses-exit: Stop refusing a failing step's commit keeps its attempts (harness.md §5, §6)", async () => {
      const s = await createSession(fixture, { model: [{ fail: { message: 'boom', transient: false } }] });
      let stopped = false;
      const port = delegate(s.port, {
        commit: async (lease, batch) => {
          if (!stopped && batch.ending === 'abandoned' && batch.end?.status === 'failed') {
            stopped = true;
            await s.client.stop({ key: 'stop-1', turnId: batch.turnId });
          }
          return await s.port.commit(lease, batch);
        },
      });
      const { launch } = await s.send('go', { run: false });
      await runWith(s, launch, port);
      const last = stepEvents(await s.events()).at(-1);
      expect(last?.end?.status).toBe('stopped');
      expect(s.model.invocations).toHaveLength(1);
      expect(last?.calls.flatMap((c) => c.attempts)).toHaveLength(1);
    });

    it('R5 extend-in-interrupted-step: an extend made in an interrupted step is stored on the turn (harness.md §7)', async () => {
      let asked = 0;
      const s = await createSession(fixture, {
        model: [{ text: 'slow', waitMs: 60_000 }, { text: 'after' }],
        configuration: { limits: { modelCalls: 0 } },
        drive: with_({
          budget: () => {
            asked += 1;
            return { extend: { modelCalls: 10 } };
          },
        }),
      });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);
      const turnId = (await s.client.execution()).turn?.id as string;
      await s.client.submit({ key: 'int-1', messages: [user('redirect')], interrupt: { turnId } });
      await s.clock.runUntil(running);
      expect(stepEvents(await s.events()).some((e) => e.ending === 'interrupted')).toBe(true);
      expect(asked).toBe(1);
    });

    it('R6 interrupted-stored-marks: an interrupted commit marks only what it sent as stored (harness.md §2, §12)', async () => {
      const later: Later = {};
      let asked = 0;
      const s = await createSession(fixture, {
        model: [{ text: 'Done.' }, { text: 'Done.' }],
        runnerLimits: { cacheBytes: 1_000_000 },
        drive: with_({
          boundary: async (_d, context) => {
            asked += 1;
            if (asked > 1) return 'proceed' as const;
            await later.s?.client.submit({
              key: 'int-1',
              messages: [user('one more thing')],
              interrupt: { turnId: context.turnId as string },
            });
            return await new Promise<never>(() => undefined);
          },
        }),
      });
      later.s = s;
      await s.send('go');
      expect(stepEvents(await s.events()).filter((e) => e.ending === 'interrupted')).toHaveLength(1);
      expect((await s.main()).at(-1)?.role).toBe('assistant');
    });

    it('R7 stale-with-expired-lease: a commit refused for an expired lease ends the hold, not a re-run (harness.md §4, §16)', async () => {
      const later: Later = {};
      const s = await createSession(fixture, {
        tools: {
          t: async (context) => {
            await later.s?.clock.sleep(35_000, context.signal);
            return { text: 't done' };
          },
        },
        model: [{ calls: [{ tool: 't', arguments: {} }] }, { text: 'final' }],
      });
      later.s = s;
      let faults = 0;
      const port = delegate(s.port, {
        renew: async (lease, input) => {
          if (faults < 3) {
            faults += 1;
            throw new Error('network down');
          }
          return await s.port.renew(lease, input);
        },
      });
      const { launch } = await s.send('go', { run: false });
      expect((await runWith(s, launch, port)).reason).toBe('lost');
      expect(s.tools.callsOf('t')).toHaveLength(1);
    });

    it('R8 interrupt-during-commit: an interrupt seen while a step commits is asked in the next step (harness.md §10)', async () => {
      const answers: string[] = [];
      const s = await createSession(fixture, {
        tools: { a: () => ({ text: 'a done' }) },
        model: [{ calls: [{ tool: 'a', arguments: {} }] }, { text: 'final' }],
        drive: with_({
          interrupt: ({ payload }) => {
            const answer = payload.stepKind === 'tools' ? 'skip' : 'abandon';
            answers.push(answer);
            return answer;
          },
        }),
      });
      let sent = false;
      const port = delegate(s.port, {
        commit: async (lease, batch) => {
          if (!sent && batch.next?.kind === 'tools') {
            sent = true;
            await s.client.submit({ key: 'int-1', messages: [user('stop that')], interrupt: { turnId: batch.turnId } });
          }
          return await s.port.commit(lease, batch);
        },
      });
      const { launch } = await s.send('go', { run: false });
      await runWith(s, launch, port);
      expect(answers).toEqual(['skip']);
      expect(s.tools.callsOf('a')).toHaveLength(0);
    });

    it('R9 withdraw-during-boundary: a steered item withdrawn while boundary is open asks again (harness.md §10)', async () => {
      let asks = 0;
      const s: Session = await createSession(fixture, {
        model: [{ text: 'one' }],
        drive: with_({
          'next-turn': (_d, context) => ({ begin: { keys: context.inbox.slice(0, 1).map((item) => item.key) } }),
          boundary: async () => {
            asks += 1;
            if (asks > 1) return 'proceed';
            await s.client.withdraw({ key: 'second' });
            await settle();
            return { steer: ['second'] };
          },
        }),
      });
      const { launch } = await s.send('first', { key: 'first', run: false });
      await s.client.submit({ key: 'second', messages: [user('second')] });
      await s.run(launch as string);
      expect((await s.turns())[0]?.status).toBe('completed');
    });

    it('R10 interrupt-of-ended-turn: an interrupt item outlives its turn as an ordinary item (harness.md §10)', async () => {
      const later: Later = {};
      let answered = 0;
      const s = await createSession(fixture, {
        tools: {
          slow: async (context) => {
            await later.s?.clock.sleep(60_000, context.signal);
            return { text: 'slow done' };
          },
        },
        model: [{ calls: [{ tool: 'slow', arguments: {} }] }, { text: 'two' }, { text: 'three' }],
        drive: with_({
          interrupt: () => {
            answered += 1;
            return 'skip';
          },
          'next-turn': (_d, context) => {
            const keys = context.inbox.filter((i) => !i.interrupt && i.key !== 'int-1').map((i) => i.key);
            return keys.length > 0 ? { begin: { keys } } : 'release';
          },
          boundary: () => 'proceed',
        }),
      });
      later.s = s;
      const { launch } = await s.send('go', { key: 'first', run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.tools.callsOf('slow').length === 1);
      const turnId = (await s.client.execution()).turn?.id as string;
      await s.client.submit({ key: 'int-1', messages: [user('redirect')], interrupt: { turnId } });
      await until(() => answered === 1);
      await s.client.submit({ key: 'second', messages: [user('second')] });
      await s.client.stop({ key: 'stop-1', turnId });
      await s.clock.runUntil(running);
      expect((await s.client.execution()).inbox.map((i) => i.key)).toEqual(['int-1']);
    });

    it('R11 interrupt-withdrawn: an interrupt withdrawn before its commit still ends the step interrupted (harness.md §12)', async () => {
      const later: Later = {};
      const s = await createSession(fixture, {
        model: [{ text: 'slow', waitMs: 60_000 }, { text: 'again' }],
        drive: with_({
          interrupt: async ({ payload }) => {
            await later.s?.client.withdraw({ key: payload.item.key });
            return 'abandon' as const;
          },
        }),
      });
      later.s = s;
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);
      const turnId = (await s.client.execution()).turn?.id as string;
      await s.client.submit({ key: 'int-1', messages: [user('redirect')], interrupt: { turnId } });
      expect((await s.clock.runUntil(running)).reason).toBe('released');
      const interrupted = stepEvents(await s.events()).find((e) => e.ending === 'interrupted');
      expect(interrupted?.consumed).toEqual([]);
      expect(interrupted?.calls.flatMap((c) => c.attempts)).toHaveLength(1);
    });

    it('R12 no-adapter: a model that cannot be resolved opens failure; the turn fails, the runner releases (harness.md §6)', async () => {
      const failures: string[] = [];
      const s = await createSession(fixture, {
        model: [{ text: 'ok' }],
        drive: with_({
          failure: ({ payload }) => {
            failures.push(payload.failure.code);
            return 'pass';
          },
        }),
      });
      const { launch } = await s.send('go', { run: false });
      const runner = s.runner({
        models: {
          resolve: async () => {
            throw new Error('no adapter for scripted-1');
          },
        },
      });
      expect((await s.run(launch as string, runner)).reason).toBe('released');
      expect(failures).toEqual(['provider']);
      expect((await s.turns())[0]?.status).toBe('failed');
    });
  });
}
