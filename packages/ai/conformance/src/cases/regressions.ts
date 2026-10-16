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
import { user } from '../testing/messages.js';
import { delegate, stepEvents, with_ } from './helpers.js';

const runWith = (s: Session, launch: LaunchKey | undefined, port: RunnerPort) =>
  s.run(launch as string, s.runner({ port }));

const scripted = { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} };

/** R1–R12: the phase-4 review's proven findings (progress.md), each kept as a case. */
export function regressionCases(fixture: ConformanceFixture): void {
  describe('Regressions', () => {
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
