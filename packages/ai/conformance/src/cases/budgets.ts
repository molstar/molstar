/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type ConformanceFixture, createSession, type Session } from '../session.js';
import { SimulatedClock, settle } from '../testing/clock.js';
import { lastStep, stepEvents, with_ } from './helpers.js';

/** Known lease and launch windows, so the cases can expire them (harness.md §4). */
const WINDOWS = { leaseMs: 30_000, launchMs: 10_000 };

/**
 * A crash: the runner stops mid-step with no commit and no release, the lease runs out, and `check`
 * launches the runner that re-runs the step (harness.md §4, §6 "lost"). Returns that launch key.
 */
async function crash(s: Session, launch: string): Promise<string> {
  const runner = s.runner();
  const running = runner.run(launch);
  await settle(40);
  runner.abort('crash');
  await s.clock.runUntil(running);
  await s.clock.advance(WINDOWS.leaseMs + 1_000);
  const checked = await s.client.check();
  expect(checked.launch).toBeDefined();
  return checked.launch as string;
}

const modelAttempts = async (s: Session) =>
  stepEvents(await s.events()).flatMap((event) =>
    event.calls.flatMap((call) => (call.kind === 'model' ? call.attempts : [])),
  );

/** B1–B2 (Budgets). */
export function budgetCases(fixture: ConformanceFixture): void {
  describe('Budgets', () => {
    // B1 (harness.md §11, Charging: every ended attempt counts once, when its step event is written)
    it('B1 usage-once: a step re-run after a crash charges only the committed attempt', async () => {
      const s = await createSession(fixture, {
        host: WINDOWS,
        model: [
          { text: 'lost', waitMs: 60_000, usage: { inputTokens: 1_000, outputTokens: 100 } },
          { text: 'kept', usage: { inputTokens: 200, outputTokens: 20 } },
        ],
      });
      const { launch } = await s.send('go', { run: false });
      await s.run(await crash(s, launch as string));
      const [turn] = await s.turns();
      expect(turn?.status).toBe('completed');
      expect(turn?.usage).toMatchObject({ modelCalls: 1, uncachedInput: 200, output: 20 });
      expect(await modelAttempts(s)).toHaveLength(1);
    });

    it('B1 missing-usage: an attempt with no reported usage is charged estimates, marked estimated', async () => {
      const s = await createSession(fixture, { model: [{ text: 'ok', usage: null }] });
      await s.send('go');
      const [turn] = await s.turns();
      expect(turn?.usage.modelCalls).toBe(1);
      expect(turn?.usage.uncachedInput).toBeGreaterThan(0);
      expect(turn?.usage.estimated).toBe(true);
      const [attempt] = await modelAttempts(s);
      expect(attempt?.usage.estimated).toBe(true);
    });

    it('B1 cache-categories: cache reads and writes are kept apart from uncached input', async () => {
      // `inputTokens` is total input, cache reads and writes included (progress.md, phase 2).
      const usage = { inputTokens: 145, outputTokens: 10, cacheReadInputTokens: 40, cacheWriteInputTokens: 5 };
      const s = await createSession(fixture, { model: [{ text: 'ok', usage }] });
      await s.send('go');
      const [turn] = await s.turns();
      expect(turn?.usage).toMatchObject({ uncachedInput: 100, cacheRead: 40, cacheWrite: 5, output: 10 });
      expect(turn?.usage.estimated).not.toBe(true);
    });

    // B2 (harness.md §11, Admission)
    it('B2 retry-after-call-limit: a retry is no new call, so it is admitted at modelCalls 1 (defect 1)', async () => {
      const s = await createSession(fixture, {
        configuration: { limits: { modelCalls: 1 } },
        model: [{ fail: { transient: true } }, { text: 'ok' }],
      });
      await s.send('go');
      const [turn] = await s.turns();
      expect(turn?.status).toBe('completed');
      expect(s.model.invocations).toHaveLength(2);
    });

    it('B2 token-shortfall: a projected input over what remains opens budget before the attempt; extend admits it', async () => {
      const asked: { dimension: string; need: number | undefined }[] = [];
      const s = await createSession(fixture, {
        configuration: { limits: { uncachedInput: 10 } },
        model: [{ text: 'ok' }],
        drive: with_({
          budget: ({ payload }) => {
            asked.push({ dimension: payload.dimension, need: payload.need });
            return { extend: { uncachedInput: 1_000_000 } };
          },
        }),
      });
      const runner = s.runner();
      await s.send('go', { runner });
      expect(asked).toHaveLength(1);
      expect(asked[0]?.dimension).toBe('uncachedInput');
      expect(asked[0]?.need).toBeGreaterThan(10);
      const order = runner.observations.map((o) =>
        o.kind === 'decision.opened' ? `opened:${o.decisionKind}` : o.kind,
      );
      expect(order.indexOf('opened:budget')).toBeGreaterThanOrEqual(0);
      expect(order.indexOf('opened:budget')).toBeLessThan(order.indexOf('attempt.started'));
      expect(s.model.invocations).toHaveLength(1);
      expect((await s.turns())[0]?.status).toBe('completed');
      expect(stepEvents(await s.events()).some((event) => event.limits?.uncachedInput === 1_000_000)).toBe(true);
    });

    it("B2 extend-persisted: an extend lost to a crash before the step's commit is asked again", async () => {
      let asks = 0;
      const s = await createSession(fixture, {
        host: WINDOWS,
        configuration: { limits: { uncachedInput: 10 } },
        model: [{ text: 'a', waitMs: 60_000 }, { text: 'b' }],
        drive: with_({
          budget: () => {
            asks += 1;
            return { extend: { uncachedInput: 1_000_000 } };
          },
        }),
      });
      const { launch } = await s.send('go', { run: false });
      const again = await crash(s, launch as string);
      expect(asks).toBe(1);
      expect((await s.client.execution()).turn?.budget.limits.uncachedInput).toBe(10);
      await s.run(again);
      expect(asks).toBe(2);
      expect((await s.turns())[0]?.status).toBe('completed');
    });

    it('B2 extend-persisted: after a commit, a new runner sees the raised limits and is not asked again', async () => {
      let asks = 0;
      let clock: Session['clock'] | undefined;
      const s = await createSession(fixture, {
        host: WINDOWS,
        configuration: { limits: { uncachedInput: 10 } },
        tools: {
          slow: async (context) => {
            await clock?.sleep(60_000, context.signal);
            return { text: 'done' };
          },
        },
        model: [{ calls: [{ tool: 'slow', arguments: {} }] }, { text: 'final' }],
        drive: with_({
          budget: () => {
            asks += 1;
            return { extend: { uncachedInput: 1_000_000 } };
          },
        }),
      });
      clock = s.clock;
      const { launch } = await s.send('go', { run: false });
      const again = await crash(s, launch as string);
      expect(asks).toBe(1);
      expect((await s.client.execution()).turn?.budget.limits.uncachedInput).toBe(1_000_000);
      await s.run(again);
      expect(asks).toBe(1);
      expect((await s.turns())[0]?.status).toBe('completed');
    });

    it('B2 end-by-mid-step: endBy passing during a tool opens budget; without an extend the step is abandoned, failed', async () => {
      const endBy = new SimulatedClock().now() + 30_000;
      let clock: Session['clock'] | undefined;
      const s = await createSession(fixture, {
        configuration: { limits: { endBy } },
        tools: {
          slow: async (context) => {
            await clock?.sleep(120_000, context.signal);
            return { text: 'late' };
          },
        },
        model: [{ calls: [{ tool: 'slow', arguments: {} }] }, { text: 'never' }],
      });
      clock = s.clock;
      const runner = s.runner();
      await s.send('go', { runner });
      const opened = runner.observations.flatMap((o) => (o.kind === 'decision.opened' ? [o.decisionKind] : []));
      expect(opened).toContain('budget');
      const last = lastStep(await s.events());
      expect(last.stepKind).toBe('tools');
      expect(last.ending).toBe('abandoned');
      expect(last.end?.status).toBe('failed');
      expect(JSON.stringify(last.end?.reason)).toContain('time');
      expect(last.calls[0]?.attempts.at(-1)?.outcome).toBe('aborted');
      expect(s.tools.callsOf('slow')[0]?.signal.aborted).toBe(true);
      expect(s.model.invocations).toHaveLength(1);
      expect((await s.turns())[0]?.status).toBe('failed');
    });
  });
}
