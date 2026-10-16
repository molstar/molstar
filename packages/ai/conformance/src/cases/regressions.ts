/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createSession, type Session } from '../session.js';
import { user } from '../testing/messages.js';
import { stepEvents, with_ } from './helpers.js';

const scripted = { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} };

/** R1–R12: the phase-4 review's proven findings (progress.md), each kept as a case. */
export function regressionCases(): void {
  describe('Regressions', () => {
    it("R3 drive-run-in-model-result: a drive's nested run does not replace the step's message (harness.md §6, §7)", async () => {
      const s = createSession({
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

    it("R4 stop-refuses-exit: Stop asked as a step fails wins over the fail, and keeps the step's attempts (harness.md §5, §6)", async () => {
      const s: Session = createSession({
        model: [{ fail: { message: 'boom', transient: false } }],
        drive: with_({
          failure: () => {
            s.harness.stop();
            return { fail: { reason: 'boom' } };
          },
        }),
      });
      await s.send('go');
      const last = stepEvents(await s.events()).at(-1);
      expect(last?.end?.status).toBe('stopped');
      expect(s.model.invocations).toHaveLength(1);
      expect(last?.calls.flatMap((c) => c.attempts)).toHaveLength(1);
    });

    it('R12 no-adapter: a model that cannot be resolved opens failure; the turn fails (harness.md §6)', async () => {
      const failures: string[] = [];
      const s = createSession({
        model: [{ text: 'ok' }],
        models: () => ({
          resolve: async () => {
            throw new Error('no adapter for scripted-1');
          },
        }),
        drive: with_({
          failure: ({ payload }) => {
            failures.push(payload.failure.code);
            return 'pass';
          },
        }),
      });
      expect((await s.send('go')).outcome).toMatchObject({ end: { status: 'failed' } });
      expect(failures).toEqual(['provider']);
      expect((await s.turns())[0]?.status).toBe('failed');
    });
  });
}
