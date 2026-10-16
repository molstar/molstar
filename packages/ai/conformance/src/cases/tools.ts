/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { AttemptFailure, Decision, StepEvent, ToolFrame } from '@molstar/ai-harness';
import { type ConformanceFixture, createSession, type Session } from '../session.js';
import { lastStep, stepEvents, with_ } from './helpers.js';

/** Handlers run inside the session they are made for; this is their way back to it. */
interface Later {
  s?: Session;
}

const nap = (later: Later, ms: number, signal: AbortSignal): Promise<void> =>
  later.s?.clock.sleep(ms, signal) ?? Promise.resolve();

/** One stage, every call at once. */
const parallel = ({ payload }: Decision<'schedule'>) => ({
  stages: [payload.calls.map((call) => call.callId)],
  maxConcurrency: payload.calls.length,
});

/** One stage per call, in call order. */
const serial = ({ payload }: Decision<'schedule'>) => ({ stages: payload.calls.map((call) => [call.callId]) });

/** Every `failure` decision's failure goes into `into`; the failure is the call's result. */
const reportFailures = (into: AttemptFailure[]) =>
  with_({
    failure: ({ payload }) => {
      into.push(payload.failure);
      return { report: {} };
    },
  });

/** Every tool result in `main`, in order. */
async function resultsOf(s: Session): Promise<{ status: string; text: string }[]> {
  return (await s.main()).flatMap((message) =>
    message.role !== 'user'
      ? []
      : message.parts.flatMap((part) =>
          part.kind === 'tool-result'
            ? [{ status: part.status, text: part.parts.map((p) => (p.kind === 'text' ? p.text : '')).join('') }]
            : [],
        ),
  );
}

const toolsStep = (events: Awaited<ReturnType<Session['events']>>): StepEvent | undefined =>
  stepEvents(events).find((event) => event.stepKind === 'tools');

/** X1–X5 (Tools). */
export function toolCases(fixture: ConformanceFixture): void {
  describe('Tools', () => {
    it('X1 denied: a denied call gives one failed result with its text and never runs (harness.md §6)', async () => {
      const s = await createSession(fixture, {
        configuration: { tools: { enabled: ['a'], preauthorized: [] } },
        tools: { a: () => ({ text: 'a done' }) },
        model: [{ calls: [{ tool: 'a', arguments: {} }] }, { text: 'final' }],
        drive: with_({ permission: () => ({ deny: 'not today' }) }),
      });
      await s.send('go');
      const results = await resultsOf(s);
      expect(results).toHaveLength(1);
      expect(results[0]?.status).toBe('failed');
      expect(results[0]?.text).toMatch(/^Denied:/);
      expect(results[0]?.text).toContain('not today');
      expect(s.tools.callsOf('a')).toHaveLength(0);
    });

    it('X1 unbound-report: a call reported at repair gives one failed result with the reported text (harness.md §8)', async () => {
      const reasons: string[] = [];
      const s = await createSession(fixture, {
        tools: { echo: () => ({ text: 'echo' }) },
        model: [{ calls: [{ tool: 'ghost', arguments: {} }] }, { text: 'final' }],
        drive: with_({
          repair: ({ payload }) => {
            reasons.push(`${payload.call.toolId}:${payload.reason}`);
            return { report: { parts: [{ kind: 'text', text: 'No tool named ghost. Use echo.' }] } };
          },
        }),
      });
      await s.send('go');
      expect(reasons).toEqual(['ghost:unknown-tool']);
      expect(await resultsOf(s)).toEqual([{ status: 'failed', text: 'No tool named ghost. Use echo.' }]);
    });

    it("X1 failure-report: a failure answered report is the call's one failed result (harness.md §6, §7)", async () => {
      const failures: AttemptFailure[] = [];
      const s = await createSession(fixture, {
        tools: { a: () => ({ status: 'failed', text: 'boom' }) },
        model: [{ calls: [{ tool: 'a', arguments: {} }] }, { text: 'final' }],
        drive: reportFailures(failures),
      });
      await s.send('go');
      expect(failures).toHaveLength(1);
      expect(s.tools.callsOf('a')).toHaveLength(1);
      expect(await resultsOf(s)).toEqual([{ status: 'failed', text: 'boom' }]);
    });

    it('X1 invalid-arguments: non-object arguments open repair and skip the handler; schema-rejected ones reach it (harness.md §8)', async () => {
      const reasons: string[] = [];
      const s = await createSession(fixture, {
        tools: {
          echo: () => ({ text: 'echo' }),
          strict: {
            definition: {
              tool: { id: 'strict' },
              name: 'strict',
              description: 'Needs a number n.',
              inputSchema: { type: 'object', properties: { n: { type: 'number' } }, required: ['n'] },
            },
            run: () => ({ text: 'normalized' }),
          },
        },
        model: [
          {
            calls: [
              { tool: 'echo', arguments: 'oops' },
              { tool: 'echo', arguments: [1, 2] },
              { tool: 'strict', arguments: {} },
            ],
          },
          { text: 'final' },
        ],
        drive: with_({
          repair: ({ payload }) => {
            reasons.push(payload.reason);
            return { report: { parts: [{ kind: 'text', text: 'arguments must be an object' }] } };
          },
        }),
      });
      await s.send('go');
      expect(reasons).toEqual(['invalid-arguments', 'invalid-arguments']);
      expect(s.tools.callsOf('echo')).toHaveLength(0);
      expect(s.tools.callsOf('strict')).toHaveLength(1);
      expect(s.tools.callsOf('strict')[0]?.arguments).toEqual({});
      expect(await resultsOf(s)).toEqual([
        { status: 'failed', text: 'arguments must be an object' },
        { status: 'failed', text: 'arguments must be an object' },
        { status: 'succeeded', text: 'normalized' },
      ]);
    });

    it('X1 skip-rest: report with a skip gives the unstarted calls failed results with that text (harness.md §7)', async () => {
      const failures: AttemptFailure[] = [];
      const s = await createSession(fixture, {
        tools: {
          a: () => ({ status: 'failed', text: 'a broke' }),
          b: () => ({ text: 'b done' }),
          c: () => ({ text: 'c done' }),
        },
        model: [
          {
            calls: [
              { tool: 'a', arguments: {} },
              { tool: 'b', arguments: {} },
              { tool: 'c', arguments: {} },
            ],
          },
          { text: 'final' },
        ],
        drive: with_({
          schedule: serial,
          failure: ({ payload }) => {
            failures.push(payload.failure);
            return { report: { skip: { text: 'Skipped: a failed' } } };
          },
        }),
      });
      await s.send('go');
      expect(failures).toHaveLength(1);
      expect(await resultsOf(s)).toEqual([
        { status: 'failed', text: 'a broke' },
        { status: 'failed', text: 'Skipped: a failed' },
        { status: 'failed', text: 'Skipped: a failed' },
      ]);
      expect(s.tools.callsOf('b')).toHaveLength(0);
      expect(s.tools.callsOf('c')).toHaveLength(0);
    });

    it('X1 no-result: a stream with no result fails the attempt stream-protocol, effect unknown (harness.md §8)', async () => {
      const failures: AttemptFailure[] = [];
      const s = await createSession(fixture, {
        tools: {
          silent: async function* (): AsyncGenerator<ToolFrame> {
            yield { kind: 'progress', value: 'working' };
          },
        },
        model: [{ calls: [{ tool: 'silent', arguments: {} }] }, { text: 'final' }],
        drive: reportFailures(failures),
      });
      await s.send('go');
      expect(failures).toHaveLength(1);
      expect(failures[0]).toMatchObject({ code: 'stream-protocol', effect: 'unknown' });
      const attempt = toolsStep(await s.events())?.calls[0]?.attempts[0];
      expect(attempt?.failure).toMatchObject({ code: 'stream-protocol', effect: 'unknown' });
      expect((await resultsOf(s))[0]?.status).toBe('failed');
    });

    it("X2 park-no-effects: with the second call's permission parked, neither handler ran (harness.md §6)", async () => {
      const s = await createSession(fixture, {
        configuration: { tools: { enabled: ['a', 'b'], preauthorized: [] } },
        tools: { a: () => ({ text: 'a done' }), b: () => ({ text: 'b done' }) },
        model: [
          {
            calls: [
              { tool: 'a', arguments: {} },
              { tool: 'b', arguments: {} },
            ],
          },
          { text: 'final' },
        ],
        drive: with_({ permission: ({ payload }) => (payload.call.toolId === 'b' ? { park: {} } : 'allow') }),
      });
      const { exit } = await s.send('go');
      expect(exit?.reason).toBe('released');
      expect((await s.client.execution()).turn?.park).toBeDefined();
      expect(lastStep(await s.events()).ending).toBe('parked');
      expect(s.tools.callsOf('a')).toHaveLength(0);
      expect(s.tools.callsOf('b')).toHaveLength(0);
    });

    it('X3 key-across-retry: attempt 2 has the same key, prior possible, attempt 2', async () => {
      let attempts = 0;
      const s = await createSession(fixture, {
        tools: {
          flaky: () => {
            attempts += 1;
            return attempts === 1 ? { status: 'failed', text: 'try again', transient: true } : { text: 'ok' };
          },
        },
        model: [{ calls: [{ tool: 'flaky', arguments: {} }] }, { text: 'final' }],
      });
      await s.send('go');
      const [first, second] = s.tools.callsOf('flaky');
      expect(first?.prior).toBe('none');
      expect(first?.attempt).toBe(1);
      expect(second?.key).toBe(first?.key);
      expect(second?.prior).toBe('possible');
      expect(second?.attempt).toBe(2);
    });

    it("X5 sibling-failure: when one call fails its siblings' signals stay unaborted (harness.md §8)", async () => {
      const later: Later = {};
      const aborted: Record<string, boolean> = {};
      const waiter = (name: string) => async (context: { signal: AbortSignal }) => {
        try {
          await nap(later, 1_000, context.signal);
        } finally {
          aborted[name] = context.signal.aborted;
        }
        return { text: `${name} done` };
      };
      const failures: AttemptFailure[] = [];
      const s = await createSession(fixture, {
        tools: { a: () => ({ status: 'failed', text: 'a failed' }), b: waiter('b'), c: waiter('c') },
        model: [
          {
            calls: [
              { tool: 'a', arguments: {} },
              { tool: 'b', arguments: {} },
              { tool: 'c', arguments: {} },
            ],
          },
          { text: 'final' },
        ],
        drive: with_({
          schedule: parallel,
          failure: ({ payload }) => {
            failures.push(payload.failure);
            return { report: {} };
          },
        }),
      });
      later.s = s;
      await s.send('go');
      expect(failures).toHaveLength(1);
      expect(aborted).toEqual({ b: false, c: false });
      expect(await resultsOf(s)).toEqual([
        { status: 'failed', text: 'a failed' },
        { status: 'succeeded', text: 'b done' },
        { status: 'succeeded', text: 'c done' },
      ]);
    });
  });
}
