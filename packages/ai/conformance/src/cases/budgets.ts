/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createSession, type Session } from '../session.js';
import { SimulatedClock } from '../testing/clock.js';
import { lastStep, stepEvents, with_ } from './helpers.js';

const modelAttempts = async (s: Session) =>
  stepEvents(await s.events()).flatMap((event) =>
    event.calls.flatMap((call) => (call.kind === 'model' ? call.attempts : [])),
  );

/** B1–B2 (Budgets). */
export function budgetCases(): void {
  describe('Budgets', () => {
    // B1 (harness.md §11, Charging: every ended attempt counts once, when its step event is written)
    it('B1 missing-usage: an attempt with no reported usage is charged estimates, marked estimated', async () => {
      const s = createSession({ model: [{ text: 'ok', usage: null }] });
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
      const s = createSession({ model: [{ text: 'ok', usage }] });
      await s.send('go');
      const [turn] = await s.turns();
      expect(turn?.usage).toMatchObject({ uncachedInput: 100, cacheRead: 40, cacheWrite: 5, output: 10 });
      expect(turn?.usage.estimated).not.toBe(true);
    });

    // B2 (harness.md §11, Admission)
    it('B2 token-shortfall: a projected input over what remains opens budget before the attempt; extend admits it', async () => {
      const asked: { dimension: string; need: number | undefined }[] = [];
      const s = createSession({
        configuration: { limits: { uncachedInput: 10 } },
        model: [{ text: 'ok' }],
        drive: with_({
          budget: ({ payload }) => {
            asked.push({ dimension: payload.dimension, need: payload.need });
            return { extend: { uncachedInput: 1_000_000 } };
          },
        }),
      });
      await s.send('go');
      expect(asked).toHaveLength(1);
      expect(asked[0]?.dimension).toBe('uncachedInput');
      expect(asked[0]?.need).toBeGreaterThan(10);
      const order = s.observations.map((o) => (o.kind === 'decision.opened' ? `opened:${o.decisionKind}` : o.kind));
      expect(order.indexOf('opened:budget')).toBeGreaterThanOrEqual(0);
      expect(order.indexOf('opened:budget')).toBeLessThan(order.indexOf('attempt.started'));
      expect(s.model.invocations).toHaveLength(1);
      expect((await s.turns())[0]?.status).toBe('completed');
      expect(stepEvents(await s.events()).some((event) => event.limits?.uncachedInput === 1_000_000)).toBe(true);
    });

    it('B2 end-by-mid-step: endBy passing during a tool opens budget; without an extend the step is abandoned, failed', async () => {
      const endBy = new SimulatedClock().now() + 30_000;
      let clock: Session['clock'] | undefined;
      const s = createSession({
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
      await s.send('go');
      const opened = s.observations.flatMap((o) => (o.kind === 'decision.opened' ? [o.decisionKind] : []));
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
