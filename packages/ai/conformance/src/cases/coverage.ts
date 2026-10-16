/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type CallView, drive, type RunResult, type RunSpec, type TurnEnd } from '@molstar/ai-harness';
import { createSession, type Session, type SessionSetup } from '../session.js';
import { textOf, toolDefinition, user } from '../testing/messages.js';
import { controls, expectWith, lastStep, stepEvents, TURN, until, with_ } from './helpers.js';

/** A nested run over the given tools, all preauthorized, on the scripted model. */
const nestedRun = (...tools: string[]): RunSpec => ({
  input: [user('work on this')],
  tools: { enabled: tools, preauthorized: tools, definitions: tools.map((id) => toolDefinition(id)) },
  model: TURN.model,
  limits: {},
});

/** A nested run whose tools need a `permission` answer. */
const guardedRun = (...tools: string[]): RunSpec => ({
  ...nestedRun(...tools),
  tools: { enabled: tools, preauthorized: [], definitions: tools.map((id) => toolDefinition(id)) },
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

/** How the first turn that ended did. */
async function endOf(s: Session): Promise<TurnEnd | undefined> {
  return stepEvents(await s.events()).find((event) => event.end !== undefined)?.end;
}

const statuses = async (s: Session): Promise<string[]> => (await s.turns()).map((turn) => turn.status);

/** A call bound to nothing: which reason `repair` was opened with, and whether the handler ran. */
async function unbound(tool: string) {
  const reasons: string[] = [];
  const s = createSession({
    declare: ['echo', 'off', 'headless'],
    deferred: ['later'],
    configuration: {
      tools: { enabled: ['echo', 'later', 'headless'], preauthorized: ['echo', 'later', 'off', 'headless'] },
    },
    tools: {
      echo: () => ({ text: 'echo' }),
      later: () => ({ text: 'later ran' }),
      off: () => ({ text: 'off ran' }),
      // A definition and no handler.
      headless: {},
    },
    model: [{ calls: [{ tool, arguments: {} }] }, { text: 'final' }],
    drive: with_({
      repair: ({ payload }) => {
        reasons.push(`${payload.call.toolId}:${payload.reason}`);
        return { report: { parts: [{ kind: 'text', text: `Cannot call ${payload.call.toolId}.` }] } };
      },
    }),
  });
  await s.send('go');
  return { s, reasons };
}

/** Coverage cases: rules of the harness (harness.md) that no other case exercises. */
export function coverageCases(): void {
  describe('Coverage', () => {
    // V1. Unbound reasons (harness.md §8).

    it('V1 not-loaded: a call to a tool not loaded opens repair with not-loaded and no handler runs (harness.md §8)', async () => {
      const { s, reasons } = await unbound('later');
      expect(reasons).toEqual(['later:not-loaded']);
      expect(s.tools.calls).toHaveLength(0);
      expect(await resultsOf(s)).toEqual([{ status: 'failed', text: 'Cannot call later.' }]);
      expect((await s.turns())[0]?.status).toBe('completed');
    });

    // V2. A drive that throws (harness.md §7).

    it('V2 drive-throws: a drive that throws ends the turn failed with drive-error; the next prompt starts a new turn (harness.md §7)', async () => {
      let asked = 0;
      const s = createSession({
        model: [{ text: 'one' }, { text: 'two' }],
        drive: with_({
          boundary: () => {
            asked += 1;
            if (asked === 1) throw new Error('the drive broke');
            return 'proceed';
          },
        }),
      });
      await s.send('go');
      expect(await statuses(s)).toEqual(['failed']);
      const end = await endOf(s);
      expect(end?.status).toBe('failed');
      expect(JSON.stringify(end?.reason)).toContain('drive-error');
      expect(lastStep(await s.events()).ending).toBe('abandoned');
      expect((await s.harness.execution()).turn).toBeUndefined();
      expect((await s.main()).map((message) => message.role)).toEqual(['user']);

      await s.send('again');
      expect(await statuses(s)).toEqual(['completed', 'failed']);
      expect((await s.main()).map(textOf)).toContain('two');
    });

    // V3. Reset (harness.md §5).

    it('V3 reset-clears-park: a parked turn, then reset: no active turn, nothing re-runs, and a later prompt runs (harness.md §5)', async () => {
      let asked = 0;
      const s = createSession({
        tools: { t: () => ({ text: 'done' }) },
        configuration: { tools: { enabled: ['t'], preauthorized: [] } },
        model: [{ calls: [{ tool: 't', arguments: {} }] }, { text: 'after the reset' }],
        drive: with_({
          permission: () => {
            asked += 1;
            return { park: { note: 'ask the person' } };
          },
        }),
      });
      const first = await s.send('go');
      expect(first.outcome).toMatchObject({ parked: { note: 'ask the person' } });
      expect((await s.harness.execution()).turn?.park).toBeDefined();

      s.harness.reset();

      const execution = await s.harness.execution();
      expect(execution.turn).toBeUndefined();
      expect(await statuses(s)).toEqual(['reset']);
      expect(controls(await s.events(), 'reset')).toHaveLength(1);
      // The parked turn is gone for good: resuming it is refused.
      expect(() => s.harness.resume()).toThrow(expect.objectContaining({ code: 'invalid' }));

      await s.send('again');
      expect(asked).toBe(1);
      expect(s.tools.calls).toHaveLength(0);
      expect(await statuses(s)).toEqual(['completed', 'reset']);
      expect((await s.main()).map(textOf)).toContain('after the reset');
    });

    // V4. Nested runs (harness.md §7, §9).

    it("V4 nested-fail: a nested drive's fail ends only that run: the parent call settles and the turn completes (harness.md §9)", async () => {
      const seen: { result?: RunResult } = {};
      const s = createSession({
        declare: ['outer'],
        tools: {
          outer: async (context) => {
            seen.result = await context.run({ ...nestedRun(), drive: drive({ boundary: () => ({ fail: 'nope' }) }) });
            return { text: seen.result.status };
          },
        },
        model: [{ calls: [{ tool: 'outer', arguments: {} }] }, { text: 'nested says hi' }, { text: 'final' }],
      });
      await s.send('go');

      expect(seen.result?.status).toBe('failed');
      expect(seen.result?.failure).toBe('nope');
      expect(await resultsOf(s)).toEqual([{ status: 'succeeded', text: 'failed' }]);
      expect(await statuses(s)).toEqual(['completed']);
      expect((await s.main()).map(textOf)).toContain('final');
      expect(s.model.remaining()).toBe(0);
    });

    it('V4 nested-stop: stop from a nested drive stops the whole turn (harness.md §9)', async () => {
      const s = createSession({
        declare: ['outer'],
        tools: {
          outer: async (context) => {
            try {
              await context.run({ ...nestedRun(), drive: drive({ boundary: () => 'stop' }) });
            } catch {
              // The step is ending; the handler is aborted.
            }
            return { text: 'outer returned' };
          },
        },
        model: [{ calls: [{ tool: 'outer', arguments: {} }] }, { text: 'nested says hi' }, { text: 'never' }],
      });
      await s.send('go');

      expect(await statuses(s)).toEqual(['stopped']);
      expect(lastStep(await s.events())).toMatchObject({ ending: 'abandoned', end: { status: 'stopped' } });
      expect((await s.harness.execution()).turn).toBeUndefined();
      expect(s.model.remaining()).toBe(1);
      // The step ended the turn: no control event did.
      expect(controls(await s.events(), 'stopped')).toHaveLength(0);
    });

    it('V4 nested-pass: a nested drive that passes lets the outer drive answer, innermost first (harness.md §7)', async () => {
      const order: string[] = [];
      const s = createSession({
        declare: ['outer'],
        tools: {
          outer: async (context) => ({
            text: (
              await context.run({
                ...guardedRun('guarded'),
                drive: drive({
                  permission: ({ payload }) => {
                    order.push(`nested:${payload.call.toolId}`);
                    return 'pass';
                  },
                }),
              })
            ).status,
          }),
          guarded: () => ({ text: 'guarded ran' }),
        },
        model: [
          { calls: [{ tool: 'outer', arguments: {} }] },
          { calls: [{ tool: 'guarded', arguments: {} }] },
          { text: 'nested done' },
          { text: 'final' },
        ],
        drive: with_({
          permission: ({ payload }) => {
            order.push(`outer:${payload.call.toolId}`);
            return payload.call.toolId === 'guarded' ? { deny: 'the app says no' } : 'allow';
          },
        }),
      });
      await s.send('go');

      expect(order).toEqual(['nested:guarded', 'outer:guarded']);
      expect(s.tools.callsOf('guarded')).toHaveLength(0);
      // The nested model saw the outer drive's denial as the call's result.
      expect(JSON.stringify(s.model.projections[2]?.messages)).toContain('Denied: the app says no');
      expect(await statuses(s)).toEqual(['completed']);
    });

    // V5. Preconditions (harness.md §12).

    it('V5 resume-not-parked: resume of a turn that is not parked is refused invalid and writes nothing (harness.md §5)', async () => {
      const s = createSession({ model: [{ text: 'slow', waitMs: 60_000 }] });
      const invalid = expect.objectContaining({ code: 'invalid' });
      // No turn at all.
      expect(() => s.harness.resume()).toThrow(invalid);
      const { handle } = await s.send('go', { wait: false });
      await until(() => s.model.invocations.length === 1);
      const before = { events: (await s.events()).length, seq: (await s.harness.execution()).seq };

      expect(() => s.harness.resume()).toThrow(invalid);

      expect((await s.events()).length).toBe(before.events);
      expect((await s.harness.execution()).seq).toBe(before.seq);
      await s.run(handle);
      // An ended turn is not parked either.
      expect(() => s.harness.resume()).toThrow(invalid);
      expect(controls(await s.events(), 'resumed')).toHaveLength(0);
      expect(await statuses(s)).toEqual(['completed']);
    });

    // V6. Budgets (harness.md §11).

    it("V6 tool-calls-limit: the toolCalls limit opens budget before the call's first attempt (harness.md §11)", async () => {
      const asked: { dimension: string; need: number | undefined }[] = [];
      const s = createSession({
        configuration: { limits: { toolCalls: 1 } },
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
        drive: with_({
          budget: ({ payload }) => {
            asked.push({ dimension: payload.dimension, need: payload.need });
            return { extend: { toolCalls: 10 } };
          },
        }),
      });
      await s.send('go');

      expect(asked).toEqual([{ dimension: 'toolCalls', need: 1 }]);
      const order = s.observations;
      const opened = order.findIndex((o) => o.kind === 'decision.opened' && o.decisionKind === 'budget');
      const startedOf = (tool: string) =>
        order.findIndex((o) => o.kind === 'attempt.started' && o.callId === s.tools.callsOf(tool)[0]?.callId);
      expect(opened).toBeGreaterThanOrEqual(0);
      expect(startedOf('a')).toBeLessThan(opened);
      expect(startedOf('b')).toBeGreaterThan(opened);
      expect(s.tools.calls.map((call) => call.toolId)).toEqual(['a', 'b']);
      expect(await statuses(s)).toEqual(['completed']);
    });

    // V7. Defaults and the schedule (harness.md §6, §7).

    it('V7 schedule-invalid: a schedule that omits, repeats, invents or adds a denied call fails invalid-answer (harness.md §7)', async () => {
      type Stages = (calls: readonly CallView[], denied: string) => string[][];
      const variants: [string, Stages][] = [
        ['omits a call', (calls) => [[calls[0]?.callId as string]]],
        [
          'repeats a call',
          (calls) => [[calls[0]?.callId as string], [calls[0]?.callId as string], [calls[1]?.callId as string]],
        ],
        ['names an unknown call', (calls) => [calls.map((call) => call.callId), ['call-ghost']]],
        ['adds a denied call', (calls, denied) => [[...calls.map((call) => call.callId), denied]]],
      ];
      for (const [name, stages] of variants) {
        const denied: { id?: string } = {};
        const setup: SessionSetup = {
          configuration: { tools: { enabled: ['a', 'b', 'c'], preauthorized: [] } },
          tools: { a: () => ({ text: 'a' }), b: () => ({ text: 'b' }), c: () => ({ text: 'c' }) },
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
            permission: ({ payload }) => {
              if (payload.call.toolId !== 'c') return 'allow';
              denied.id = payload.call.callId;
              return { deny: 'not c' };
            },
            schedule: ({ payload }) => ({ stages: stages(payload.calls, denied.id as string) }),
          }),
        };
        const s = createSession(setup);
        await s.send('go');
        expectWith(await statuses(s), name).toEqual(['failed']);
        expectWith(JSON.stringify((await endOf(s))?.reason), name).toContain('invalid-answer');
        expectWith(s.tools.calls, name).toHaveLength(0);
      }
    });

    it('V7 default-model-retry: defaultDrive retries a transient model failure three attempts with growing waits, then fails the turn (harness.md §7)', async () => {
      const s = createSession({
        model: [
          { fail: { transient: true } },
          { fail: { transient: true } },
          { fail: { transient: true } },
          { text: 'never' },
        ],
      });
      await s.send('go');

      expect(s.model.invocations).toHaveLength(3);
      expect(s.model.remaining()).toBe(1);
      expect(await statuses(s)).toEqual(['failed']);
      const call = stepEvents(await s.events())
        .flatMap((event) => event.calls)
        .find((c) => c.kind === 'model');
      const attempts = call?.attempts ?? [];
      expect(attempts.map((a) => a.outcome)).toEqual(['failed', 'failed', 'failed']);
      const gaps = [1, 2].map((i) => (attempts[i]?.startedAt ?? 0) - (attempts[i - 1]?.endedAt ?? 0));
      expect(gaps[0]).toBeGreaterThan(0);
      expect(gaps[1]).toBeGreaterThan(gaps[0] as number);
      expect(lastStep(await s.events())).toMatchObject({ ending: 'abandoned', end: { status: 'failed' } });
    });

    it('V7 default-tool-retry: defaultDrive retries a transient failed result three attempts, then reports it; a plain failure is reported at once (harness.md §7)', async () => {
      const s = createSession({
        tools: {
          flaky: () => ({ status: 'failed', text: 'flaky', transient: true }),
          broken: () => ({ status: 'failed', text: 'broken' }),
        },
        model: [
          {
            calls: [
              { tool: 'flaky', arguments: {} },
              { tool: 'broken', arguments: {} },
            ],
          },
          { text: 'final' },
        ],
      });
      await s.send('go');

      expect(s.tools.callsOf('flaky').map((call) => [call.attempt, call.prior])).toEqual([
        [1, 'none'],
        [2, 'possible'],
        [3, 'possible'],
      ]);
      expect(s.tools.callsOf('broken')).toHaveLength(1);
      expect(await resultsOf(s)).toEqual([
        { status: 'failed', text: 'flaky' },
        { status: 'failed', text: 'broken' },
      ]);
      expect(await statuses(s)).toEqual(['completed']);
    });

    // V8. The model result (harness.md §6).

    it('V8 filtered: a filtered finish opens model-result too, and accept completes the turn (harness.md §6)', async () => {
      const seen: string[] = [];
      const s = createSession({
        model: [{ text: 'partial', finish: 'filtered' }],
        drive: with_({
          'model-result': ({ payload }) => {
            seen.push(payload.finish);
            return 'accept';
          },
        }),
      });
      await s.send('go');

      expect(seen).toEqual(['filtered']);
      expect((await s.main()).map(textOf)).toEqual(['go', 'partial']);
      expect(await statuses(s)).toEqual(['completed']);
    });
  });
}
