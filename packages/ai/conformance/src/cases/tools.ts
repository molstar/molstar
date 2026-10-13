/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message } from '@molstar/ai-agent';
import type {
  AttemptFailure,
  Decision,
  DriveHandlers,
  LaunchKey,
  RunnerExit,
  RunSpec,
  StepEvent,
  ToolFrame,
} from '@molstar/ai-harness';
import { type ConformanceFixture, createSession, type Session } from '../session.js';
import { settle } from '../testing/clock.js';
import { toolDefinition, user } from '../testing/messages.js';
import { lastStep, stepEvents, with_ } from './helpers.js';

/** Handlers run inside the session they are made for; this is their way back to it. */
interface Later {
  s?: Session;
}

const nap = (later: Later, ms: number, signal: AbortSignal): Promise<void> =>
  later.s?.clock.sleep(ms, signal) ?? Promise.resolve();

/** Lets real work proceed until `condition` holds. */
async function until(condition: () => boolean, limit = 500): Promise<void> {
  for (let i = 0; i < limit && !condition(); i += 1) await settle();
  if (!condition()) throw new Error('until: the condition never held');
}

/** A crash: the runner stops where it is and writes nothing; its lease is left to run out. */
async function crash(s: Session, launch: LaunchKey, reached: () => boolean): Promise<void> {
  const runner = s.runner();
  const running = runner.run(launch);
  await until(reached);
  runner.abort('crash');
  await s.clock.runUntil(running);
}

/** The lease runs out, `check` launches, and a new runner takes the turn over. */
async function recover(s: Session): Promise<RunnerExit> {
  await s.clock.advance(120_000);
  const { launch } = await s.client.check();
  if (launch === undefined) throw new Error('no launch after the lease expired');
  return await s.run(launch);
}

/** A nested run over the given tools, all preauthorized, on the scripted model. */
const nestedRun = (...tools: string[]): RunSpec => ({
  input: [user('work on this')],
  tools: { enabled: tools, preauthorized: tools, definitions: tools.map((id) => toolDefinition(id)) },
  model: { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} },
  limits: {},
});

/** One stage, every call at once. */
const parallel = ({ payload }: Decision<'schedule'>) => ({
  stages: [payload.calls.map((call) => call.callId)],
  maxConcurrency: payload.calls.length,
});

/** One stage per call, in call order. */
const serial = ({ payload }: Decision<'schedule'>) => ({ stages: payload.calls.map((call) => [call.callId]) });

/** Steers the interrupting item (key `int`) while it waits. */
const steerInterrupt: Required<DriveHandlers>['boundary'] = (_decision, context) =>
  context.inbox.some((item) => item.key === 'int') ? { steer: ['int'] } : 'proceed';

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

    it('X1 skipped: an interrupt answered skip gives each unstarted call one failed result (harness.md §6, §10)', async () => {
      const later: Later = {};
      const s = await createSession(fixture, {
        tools: {
          a: async (context) => {
            await later.s?.client.submit({
              key: 'int',
              messages: [user('hold on')],
              interrupt: { turnId: context.turnId },
            });
            await nap(later, 1_000, context.signal);
            return { text: 'a done' };
          },
          b: () => ({ text: 'b done' }),
        },
        model: [
          {
            calls: [
              { tool: 'a', arguments: {} },
              { tool: 'b', arguments: {} },
            ],
          },
          { text: 'after' },
        ],
        drive: with_({ schedule: serial, interrupt: () => 'skip', boundary: steerInterrupt }),
      });
      later.s = s;
      await s.send('go');
      const results = await resultsOf(s);
      expect(results).toHaveLength(2);
      expect(results[1]?.status).toBe('failed');
      expect(results[1]?.text).toContain('Skipped: interrupted');
      expect(s.tools.callsOf('b')).toHaveLength(0);
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

    it("X1 attempt-times: each attempt records its own start and end, not the settlement's (harness.md §3)", async () => {
      const later: Later = {};
      const s = await createSession(fixture, {
        tools: {
          a: async (context) => {
            await nap(later, 100, context.signal);
            return { text: 'a' };
          },
          b: async (context) => {
            await nap(later, 300, context.signal);
            return { text: 'b' };
          },
        },
        model: [
          {
            calls: [
              { tool: 'a', arguments: {} },
              { tool: 'b', arguments: {} },
            ],
          },
          { text: 'final' },
        ],
        drive: with_({ schedule: parallel }),
      });
      later.s = s;
      await s.send('go');
      const step = toolsStep(await s.events());
      const attemptOf = (tool: string) => {
        const callId = s.tools.callsOf(tool)[0]?.callId;
        const attempt = step?.calls.find((call) => call.id === callId)?.attempts[0];
        if (attempt === undefined) throw new Error(`no attempt recorded for ${tool}`);
        return attempt;
      };
      const [a, b] = [attemptOf('a'), attemptOf('b')];
      expect(a.endedAt - a.startedAt).toBeGreaterThanOrEqual(100);
      expect(a.endedAt - a.startedAt).toBeLessThan(300);
      expect(b.endedAt - b.startedAt).toBeGreaterThanOrEqual(300);
      expect(a.endedAt).toBeLessThan(b.endedAt);
      expect(a.endedAt).toBeLessThan(step?.at ?? 0);
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

    it('X1 skip-calls: the payload lists the unstarted calls; a skip naming some skips only those (harness.md §7)', async () => {
      const unstarted: string[][] = [];
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
          // The drive's view of what depends on a: c (b does not).
          failure: ({ payload }) => {
            const rest = payload.unstarted ?? [];
            unstarted.push(rest.map((call) => call.toolId));
            const c = rest.find((call) => call.toolId === 'c');
            return { report: { skip: { text: 'Skipped: it needed a', calls: c ? [c.callId] : [] } } };
          },
        }),
      });
      await s.send('go');
      expect(unstarted).toEqual([['b', 'c']]);
      expect(await resultsOf(s)).toEqual([
        { status: 'failed', text: 'a broke' },
        { status: 'succeeded', text: 'b done' },
        { status: 'failed', text: 'Skipped: it needed a' },
      ]);
      expect(s.tools.callsOf('b')).toHaveLength(1);
      expect(s.tools.callsOf('c')).toHaveLength(0);
    });

    it('X1 skip-calls-invalid: a skip naming a call that is not unstarted fails invalid-answer (harness.md §7)', async () => {
      const s = await createSession(fixture, {
        tools: { a: () => ({ status: 'failed', text: 'a broke' }), b: () => ({ text: 'b done' }) },
        model: [
          {
            calls: [
              { tool: 'a', arguments: {} },
              { tool: 'b', arguments: {} },
            ],
          },
          { text: 'final' },
        ],
        drive: with_({
          schedule: serial,
          // The failed call itself has started: it cannot be skipped.
          failure: ({ payload }) => ({ report: { skip: { text: 'no', calls: [payload.call?.callId as string] } } }),
        }),
      });
      await s.send('go');
      expect(JSON.stringify(lastStep(await s.events()).end?.reason)).toContain('invalid-answer');
      expect(s.tools.callsOf('b')).toHaveLength(0);
    });

    it('X1 extra-frame: a stream with two results fails the attempt stream-protocol, effect unknown (harness.md §8)', async () => {
      const failures: AttemptFailure[] = [];
      const s = await createSession(fixture, {
        tools: {
          twice: async function* (): AsyncGenerator<ToolFrame> {
            yield { kind: 'result', status: 'succeeded', parts: [{ kind: 'text', text: 'one' }] };
            yield { kind: 'result', status: 'succeeded', parts: [{ kind: 'text', text: 'two' }] };
          },
        },
        model: [{ calls: [{ tool: 'twice', arguments: {} }] }, { text: 'final' }],
        drive: reportFailures(failures),
      });
      await s.send('go');
      expect(failures).toHaveLength(1);
      expect(failures[0]).toMatchObject({ code: 'stream-protocol', effect: 'unknown' });
      const attempt = toolsStep(await s.events())?.calls[0]?.attempts[0];
      expect(attempt?.outcome).toBe('failed');
      expect(attempt?.failure).toMatchObject({ code: 'stream-protocol', effect: 'unknown' });
      expect((await resultsOf(s))[0]?.status).toBe('failed');
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

    it('X3 key-across-crash: after a crash mid-call the re-run gets the same key and prior possible (harness.md §8)', async () => {
      const later: Later = {};
      let runs = 0;
      let reached = false;
      const s = await createSession(fixture, {
        host: { leaseMs: 30_000 },
        tools: {
          work: async (context) => {
            runs += 1;
            if (runs === 1) {
              reached = true;
              await nap(later, 600_000, context.signal);
            }
            return { text: 'done' };
          },
        },
        model: [{ calls: [{ tool: 'work', arguments: {} }] }, { text: 'final' }],
      });
      later.s = s;
      const { launch } = await s.send('go', { run: false });
      await crash(s, launch as string, () => reached);
      await recover(s);
      const [first, second] = s.tools.callsOf('work');
      expect(s.tools.callsOf('work')).toHaveLength(2);
      expect(first?.prior).toBe('none');
      expect(second?.key).toBe(first?.key);
      expect(second?.prior).toBe('possible');
      expect((await s.turns())[0]?.status).toBe('completed');
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

    it('X3 nested-key-rerun: a nested call gets the same key when its parent step re-runs (harness.md §8)', async () => {
      const later: Later = {};
      let runs = 0;
      let reached = false;
      const s = await createSession(fixture, {
        host: { leaseMs: 30_000 },
        declare: ['outer'],
        tools: {
          outer: async (context) => {
            runs += 1;
            const result = await context.run(nestedRun('inner'));
            if (runs === 1) {
              reached = true;
              await nap(later, 600_000, context.signal);
            }
            return { text: result.status };
          },
          inner: () => ({ text: 'inner done' }),
        },
        model: [
          { calls: [{ tool: 'outer', arguments: {} }] },
          { calls: [{ tool: 'inner', arguments: { q: 1 } }] },
          { text: 'nested done' },
          { calls: [{ tool: 'inner', arguments: { q: 1 } }] },
          { text: 'nested done' },
          { text: 'final' },
        ],
      });
      later.s = s;
      const { launch } = await s.send('go', { run: false });
      await crash(s, launch as string, () => reached);
      await recover(s);
      const [first, second] = s.tools.callsOf('inner');
      expect(s.tools.callsOf('inner')).toHaveLength(2);
      expect(first?.path).toHaveLength(1);
      expect(second?.key).toBe(first?.key);
      expect(second?.callId).not.toBe(first?.callId);
      expect(second?.prior).toBe('possible');
    });

    it('X3 nested-identical-calls: two identical calls in one nested run get different keys (harness.md §8)', async () => {
      const s = await createSession(fixture, {
        declare: ['outer'],
        tools: {
          outer: async (context) => ({ text: (await context.run(nestedRun('inner'))).status }),
          inner: () => ({ text: 'inner done' }),
        },
        model: [
          { calls: [{ tool: 'outer', arguments: {} }] },
          {
            calls: [
              { tool: 'inner', arguments: { q: 1 } },
              { tool: 'inner', arguments: { q: 1 } },
            ],
          },
          { text: 'nested done' },
          { text: 'final' },
        ],
      });
      await s.send('go');
      const [first, second] = s.tools.callsOf('inner');
      expect(s.tools.callsOf('inner')).toHaveLength(2);
      expect(first?.arguments).toEqual(second?.arguments);
      expect(second?.key).not.toBe(first?.key);
    });

    it('X4 orphan-result: a contribution with a tool-result and no call fails invalid-contribution (harness.md §8)', async () => {
      const orphan: Message = {
        role: 'user',
        parts: [
          {
            kind: 'tool-result',
            requestId: 'ghost',
            tool: { id: 'ghost' },
            status: 'succeeded',
            parts: [{ kind: 'text', text: 'from nowhere' }],
          },
        ],
      };
      const failures: AttemptFailure[] = [];
      const s = await createSession(fixture, {
        tools: { a: () => ({ text: 'a done', contributions: [orphan] }) },
        model: [{ calls: [{ tool: 'a', arguments: {} }] }, { text: 'final' }],
        drive: reportFailures(failures),
      });
      await s.send('go');
      expect(failures).toHaveLength(1);
      expect(failures[0]?.code).toBe('invalid-contribution');
      const main = await s.main();
      expect(main.some((m) => m.parts.some((p) => p.kind === 'tool-result' && p.requestId === 'ghost'))).toBe(false);
      expect((await resultsOf(s)).map((r) => r.status)).toEqual(['failed']);
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

    it('X5 interrupt-abandon: abandon aborts every in-flight call (harness.md §10)', async () => {
      const later: Later = {};
      const aborted: Record<string, boolean> = {};
      const hang = (name: string) => async (context: { signal: AbortSignal }) => {
        try {
          await nap(later, 60_000, context.signal);
        } finally {
          aborted[name] = context.signal.aborted;
        }
        return { text: `${name} done` };
      };
      const s = await createSession(fixture, {
        tools: { a: hang('a'), b: hang('b') },
        model: [
          {
            calls: [
              { tool: 'a', arguments: {} },
              { tool: 'b', arguments: {} },
            ],
          },
          { text: 'after' },
        ],
        drive: with_({ schedule: parallel, interrupt: () => 'abandon', boundary: steerInterrupt }),
      });
      later.s = s;
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.tools.callsOf('a').length + s.tools.callsOf('b').length === 2);
      const turnId = (await s.client.execution()).turn?.id as string;
      await s.client.submit({ key: 'int', messages: [user('stop that')], interrupt: { turnId } });
      await s.clock.runUntil(running);
      expect(aborted).toEqual({ a: true, b: true });
      const interrupted = stepEvents(await s.events()).find((event) => event.ending === 'interrupted');
      expect(interrupted).toBeDefined();
      const attempts = interrupted?.calls.flatMap((call) => call.attempts) ?? [];
      expect(attempts).toHaveLength(2);
      expect(attempts.map((attempt) => attempt.outcome)).toEqual(['aborted', 'aborted']);
    });

    it("X5 exit-aborts-siblings: a parked sub-call permission aborts a sibling's nested run and records its attempts (harness.md §6)", async () => {
      const later: Later = {};
      const s = await createSession(fixture, {
        declare: ['outer1', 'outer2'],
        tools: {
          outer1: async (context) => {
            const run: RunSpec = {
              ...nestedRun(),
              tools: { enabled: ['guarded'], preauthorized: [], definitions: [toolDefinition('guarded')] },
            };
            return { text: (await context.run(run)).status };
          },
          outer2: async (context) => {
            // Starts its nested run after outer1's, so the script order is fixed.
            await until(() => (later.s?.model.invocations.length ?? 0) >= 2);
            return { text: (await context.run(nestedRun())).status };
          },
          guarded: () => ({ text: 'never' }),
        },
        model: [
          {
            calls: [
              { tool: 'outer1', arguments: {} },
              { tool: 'outer2', arguments: {} },
            ],
          },
          { calls: [{ tool: 'guarded', arguments: {} }] },
          { text: 'never', waitMs: 600_000 },
        ],
        drive: with_({
          schedule: parallel,
          permission: async ({ payload }) => {
            if (payload.call.toolId !== 'guarded') return 'allow';
            // Parks only once the sibling's nested model call is in flight.
            await until(() => (later.s?.model.invocations.length ?? 0) >= 3);
            return { park: {} };
          },
        }),
      });
      later.s = s;
      const { exit } = await s.send('go');
      expect(exit?.reason).toBe('released');
      const parked = lastStep(await s.events());
      expect(parked.ending).toBe('parked');
      const sibling = s.tools.callsOf('outer2')[0];
      expect(sibling?.signal.aborted).toBe(true);
      const nested = parked.calls.find((call) => call.parentCallId === sibling?.callId);
      expect(nested?.kind).toBe('model');
      expect(nested?.attempts).toHaveLength(1);
      expect(nested?.attempts[0]?.outcome).toBe('aborted');
      expect(s.tools.callsOf('guarded')).toHaveLength(0);
    });
  });
}
