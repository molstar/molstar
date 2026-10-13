/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { digest, hashMessage, type Json, type Message } from '@molstar/ai-agent';
import {
  type CallRecord,
  type CallView,
  drive,
  type Item,
  type LaunchKey,
  MAIN,
  type Observation,
  type RunResult,
  type RunSpec,
  SCHEMA_VERSION,
  type ToolFrame,
  type TurnEnd,
  type Usage,
} from '@molstar/ai-harness';
import { describe, expect, it } from '@jest/globals';
import { type ConformanceFixture, createSession, type Session, type SessionSetup } from '../session.js';
import { settle } from '../testing/clock.js';
import { textOf, toolDefinition, user } from '../testing/messages.js';
import { controls, expectWith, lastStep, stepEvents, TURN, with_ } from './helpers.js';

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

/** The messages of a pointer's current head. */
async function messagesOf(s: Session, scope: string): Promise<Message[]> {
  const [pointer] = await s.client.pointers([scope]);
  if (pointer === undefined) return [];
  const [head] = await s.client.heads([pointer.head]);
  return head === undefined ? [] : await s.client.entries(head.entries);
}

/** How the first turn that ended did. */
async function endOf(s: Session): Promise<TurnEnd | undefined> {
  return stepEvents(await s.events()).find((event) => event.end !== undefined)?.end;
}

const statuses = async (s: Session): Promise<string[]> => (await s.turns()).map((turn) => turn.status);

const idsOf = (calls: readonly CallRecord[]): string[] => calls.map((call) => call.id);

/** A call bound to nothing: which reason `repair` was opened with, and whether the handler ran. */
async function unbound(fixture: ConformanceFixture, tool: string) {
  const reasons: string[] = [];
  const s = await createSession(fixture, {
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

/** How many calls of one stage ran at once, and in what order they started. */
async function concurrency(fixture: ConformanceFixture, maxConcurrency: number | undefined) {
  const later: Later = {};
  let running = 0;
  let peak = 0;
  const work = async (context: { signal: AbortSignal }) => {
    running += 1;
    peak = Math.max(peak, running);
    try {
      await nap(later, 100, context.signal);
    } finally {
      running -= 1;
    }
    return { text: 'ok' };
  };
  const s = await createSession(fixture, {
    tools: { a: work, b: work, c: work },
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
      schedule: ({ payload }) => ({
        stages: [payload.calls.map((call) => call.callId)],
        ...(maxConcurrency !== undefined ? { maxConcurrency } : {}),
      }),
    }),
  });
  later.s = s;
  await s.send('go');
  expect((await s.turns())[0]?.status).toBe('completed');
  return { peak, order: s.tools.calls.map((call) => call.toolId) };
}

/** Coverage cases: rules of the harness (harness.md) that no other case exercises. */
export function coverageCases(fixture: ConformanceFixture): void {
  describe('Coverage', () => {
    // V1. Unbound reasons (harness.md §8).

    for (const [tool, reason] of [
      ['later', 'not-loaded'],
      ['off', 'not-enabled'],
      ['headless', 'no-handler'],
    ] as const) {
      it(`V1 ${reason}: a call to ${tool} opens repair with ${reason} and no handler runs (harness.md §8)`, async () => {
        const { s, reasons } = await unbound(fixture, tool);
        expect(reasons).toEqual([`${tool}:${reason}`]);
        expect(s.tools.calls).toHaveLength(0);
        expect(await resultsOf(s)).toEqual([{ status: 'failed', text: `Cannot call ${tool}.` }]);
        expect((await s.turns())[0]?.status).toBe('completed');
      });
    }

    // V2. A drive that throws (harness.md §7).

    it('V2 drive-throws: a drive that throws ends the turn failed with drive-error; the next prompt starts a new turn (harness.md §7)', async () => {
      let asked = 0;
      const s = await createSession(fixture, {
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
      expect((await s.client.execution()).turn).toBeUndefined();
      expect((await s.main()).map((message) => message.role)).toEqual(['user']);

      await s.send('again');
      expect(await statuses(s)).toEqual(['completed', 'failed']);
      expect((await s.main()).map(textOf)).toContain('two');
    });

    it('V2 drive-rejects: a drive that rejects in a tools step abandons it before any call, and the next prompt runs (harness.md §7)', async () => {
      const s = await createSession(fixture, {
        tools: { echo: () => ({ text: 'echo' }) },
        model: [{ calls: [{ tool: 'echo', arguments: {} }] }, { text: 'ok' }],
        drive: with_({
          schedule: async () => {
            await settle(2);
            throw new Error('the drive rejected');
          },
        }),
      });
      await s.send('go');
      expect(await statuses(s)).toEqual(['failed']);
      const last = lastStep(await s.events());
      expect(last).toMatchObject({ stepKind: 'tools', ending: 'abandoned', end: { status: 'failed' } });
      expect(JSON.stringify(last.end?.reason)).toContain('drive-error');
      expect(s.tools.calls).toHaveLength(0);

      await s.send('again');
      expect(await statuses(s)).toEqual(['completed', 'failed']);
    });

    // V3. Reset (harness.md §5).

    it('V3 reset-clears-park: a parked turn, then reset: no active turn, nothing re-runs, and a later prompt runs (harness.md §5)', async () => {
      let asked = 0;
      const s = await createSession(fixture, {
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
      expect(first.exit?.reason).toBe('released');
      const turnId = (await s.client.execution()).turn?.id as string;
      expect((await s.client.execution()).turn?.park).toBeDefined();

      const reset = await s.client.reset({ key: 'reset-1' });

      expect(reset.launch).toBeUndefined();
      const execution = await s.client.execution();
      expect(execution.turn).toBeUndefined();
      expect(execution.runner).toEqual({ phase: 'none' });
      expect(await statuses(s)).toEqual(['reset']);
      expect(controls(await s.events(), 'reset')).toHaveLength(1);
      // The parked turn is gone for good: resuming it is refused.
      await expect(s.client.resume({ key: 'resume-1', turnId })).rejects.toMatchObject({ code: 'invalid' });

      await s.send('again');
      expect(asked).toBe(1);
      expect(s.tools.calls).toHaveLength(0);
      expect(await statuses(s)).toEqual(['completed', 'reset']);
      expect((await s.main()).map(textOf)).toContain('after the reset');
    });

    it('V3 reset-ends-blocked: a blocked turn is ended by reset, and a later prompt runs (harness.md §5)', async () => {
      const s = await createSession(fixture, {
        host: { leaseMs: 30_000 },
        model: [{ text: 'slow', waitMs: 600_000 }, { text: 'ok' }],
      });
      const { launch } = await s.send('go', { run: false });
      await crash(s, launch as string, () => s.model.invocations.length === 1);
      await s.clock.advance(120_000);
      const { launch: relaunch } = await s.client.check();
      expect(relaunch).toBeDefined();
      await expect(
        s.port.claim({ launch: relaunch as string, runnerId: 'r-other', schemaVersion: SCHEMA_VERSION + 1 }),
      ).rejects.toMatchObject({ code: 'stale-schema' });
      const blocked = (await s.client.execution()).turn;
      expect(blocked?.blocked).toBe('stale-schema');

      await s.client.reset({ key: 'reset-1' });

      const execution = await s.client.execution();
      expect(execution.turn).toBeUndefined();
      expect(execution.runner).toEqual({ phase: 'none' });
      expect(await statuses(s)).toEqual(['reset']);
      expect(controls(await s.events(), 'reset')).toHaveLength(1);

      await s.send('again');
      expect(await statuses(s)).toEqual(['completed', 'reset']);
    });

    it('V3 reset-other-pointer: reset with a head for another pointer moves it and nothing else; a second reset undoes it (harness.md §5)', async () => {
      let notes = 0;
      const s = await createSession(fixture, {
        tools: {
          note: (context) => {
            notes += 1;
            context.journal.append('other', [user(`note ${notes}`)]);
            return { text: 'noted' };
          },
        },
        model: [
          { calls: [{ tool: 'note', arguments: {} }] },
          { text: 'a' },
          { calls: [{ tool: 'note', arguments: {} }] },
          { text: 'b' },
        ],
      });
      await s.send('one');
      const [first] = await s.client.pointers(['other']);
      await s.send('two');
      const [second] = await s.client.pointers(['other']);
      expect(first?.head).toBeDefined();
      expect(second?.head).not.toBe(first?.head);
      expect((await messagesOf(s, 'other')).map(textOf)).toEqual(['note 1', 'note 2']);
      const before = await s.client.execution();
      const turns = await statuses(s);

      await s.client.reset({ key: 'reset-1', pointers: { other: first?.head as string } });

      expect((await s.client.pointers(['other']))[0]?.head).toBe(first?.head);
      expect((await messagesOf(s, 'other')).map(textOf)).toEqual(['note 1']);
      const after = await s.client.execution();
      expect(after.main).toBe(before.main);
      expect(await statuses(s)).toEqual(turns);

      // It deleted nothing: the head it left is still there to go back to.
      await s.client.reset({ key: 'reset-2', pointers: { other: second?.head as string } });
      expect((await messagesOf(s, 'other')).map(textOf)).toEqual(['note 1', 'note 2']);
    });

    // V4. Nested runs (harness.md §7, §9).

    it("V4 nested-limits: a nested run's own limits fail it while the parent call settles (harness.md §9, §11)", async () => {
      const seen: { result?: RunResult } = {};
      const s = await createSession(fixture, {
        declare: ['outer'],
        tools: {
          outer: async (context) => {
            seen.result = await context.run({ ...nestedRun('inner'), limits: { modelCalls: 1 } });
            return { text: seen.result.status };
          },
          inner: () => ({ text: 'inner done' }),
        },
        model: [
          { calls: [{ tool: 'outer', arguments: {} }] },
          { calls: [{ tool: 'inner', arguments: {} }] },
          { text: 'final' },
        ],
      });
      const observer = s.runner();
      await s.send('go', { runner: observer });

      expect(seen.result?.status).toBe('failed');
      expect(JSON.stringify(seen.result?.failure)).toContain('budget');
      const outerId = s.tools.callsOf('outer')[0]?.callId;
      const opened = observer.observations.flatMap((o) =>
        o.kind === 'decision.opened' && o.decisionKind === 'budget' ? [o.path] : [],
      );
      expect(opened).toEqual([[outerId]]);
      // The second nested model call was never made; the parent's own model went on.
      expect(s.tools.callsOf('inner')).toHaveLength(1);
      expect(s.model.invocations).toHaveLength(3);
      expect(await resultsOf(s)).toEqual([{ status: 'succeeded', text: 'failed' }]);
      expect(await statuses(s)).toEqual(['completed']);
    });

    it("V4 nested-fail: a nested drive's fail ends only that run: the parent call settles and the turn completes (harness.md §9)", async () => {
      const seen: { result?: RunResult } = {};
      const s = await createSession(fixture, {
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
      const s = await createSession(fixture, {
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
      expect((await s.client.execution()).turn).toBeUndefined();
      expect(s.model.remaining()).toBe(1);
      expect(controls(await s.events(), 'stop-requested')).toHaveLength(0);
    });

    it('V4 nested-pass: a nested drive that passes lets the outer drive answer, innermost first (harness.md §7)', async () => {
      const order: string[] = [];
      const s = await createSession(fixture, {
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

    it('V4 nested-permission: a nested permission drive auto-approves a sub-call the outer drive would deny (harness.md §7)', async () => {
      const asked: string[] = [];
      const s = await createSession(fixture, {
        declare: ['outer'],
        tools: {
          outer: async (context) => ({
            text: (await context.run({ ...guardedRun('guarded'), drive: drive({ permission: () => 'allow' }) })).status,
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
            asked.push(payload.call.toolId);
            return { deny: 'the app says no' };
          },
        }),
      });
      await s.send('go');

      expect(s.tools.callsOf('guarded')).toHaveLength(1);
      expect(asked).toEqual([]);
      expect(await resultsOf(s)).toEqual([{ status: 'succeeded', text: 'completed' }]);
      expect(await statuses(s)).toEqual(['completed']);
    });

    // V5. Preconditions (harness.md §12).

    it('V5 resume-not-parked: resume of a turn that is not parked is refused invalid and writes nothing (harness.md §12)', async () => {
      const s = await createSession(fixture, { model: [{ text: 'slow', waitMs: 60_000 }] });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);
      const turnId = (await s.client.execution()).turn?.id as string;
      const before = { events: (await s.events()).length, seq: (await s.client.execution()).seq };

      await expect(s.client.resume({ key: 'resume-1', turnId })).rejects.toMatchObject({ code: 'invalid' });
      await expect(s.client.resume({ key: 'resume-2', turnId: 'turn-ghost' })).rejects.toMatchObject({
        code: 'invalid',
      });

      expect((await s.events()).length).toBe(before.events);
      expect((await s.client.execution()).seq).toBe(before.seq);
      await s.clock.runUntil(running);
      // An ended turn is not parked either.
      await expect(s.client.resume({ key: 'resume-3', turnId })).rejects.toMatchObject({ code: 'invalid' });
      expect(controls(await s.events(), 'resumed')).toHaveLength(0);
      expect(await statuses(s)).toEqual(['completed']);
    });

    it('V5 withdraw-missing: withdraw of a missing item returns false and writes no event (harness.md §10, §12)', async () => {
      const s = await createSession(fixture);
      await s.send('go', { key: 'it1', run: false });
      const before = { events: await s.events(), seq: (await s.client.execution()).seq };

      const missing = await s.client.withdraw({ key: 'ghost' });

      expect(missing.withdrawn).toBe(false);
      expect(await s.events()).toEqual(before.events);
      expect((await s.client.execution()).seq).toBe(before.seq);

      expect((await s.client.withdraw({ key: 'it1' })).withdrawn).toBe(true);
      const after = await s.events();
      expect(controls(after, 'withdrawn')).toHaveLength(1);
      // A repeat is a missing item too.
      expect((await s.client.withdraw({ key: 'it1' })).withdrawn).toBe(false);
      expect(await s.events()).toEqual(after);
    });

    it('V5 putblobs-expired: putBlobs with an expired lease is refused stale and stores nothing (harness.md §12)', async () => {
      const s = await createSession(fixture);
      const { launch } = await s.client.submit({ key: 'item-1', messages: [user('hi')] });
      const { lease } = await s.port.claim({
        launch: launch as string,
        runnerId: 'r1',
        schemaVersion: SCHEMA_VERSION,
      });
      const live: Json = { say: 'live' };
      const liveRef = await digest(live);
      await s.port.putBlobs(lease, [{ ref: liveRef, body: live }]);
      expect(await s.client.blobs([liveRef])).toEqual([live]);

      await s.clock.advance(lease.until - s.clock.now() + 1);
      const late: Json = { say: 'late' };
      const lateRef = await digest(late);
      await expect(s.port.putBlobs(lease, [{ ref: lateRef, body: late }])).rejects.toMatchObject({ code: 'stale' });
      expect(await s.client.blobs([liveRef])).toEqual([live]);
    });

    // V6. Budgets (harness.md §11).

    it("V6 tool-calls-limit: the toolCalls limit opens budget before the call's first attempt (harness.md §11)", async () => {
      const asked: { dimension: string; need: number | undefined }[] = [];
      const s = await createSession(fixture, {
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
      const runner = s.runner();
      await s.send('go', { runner });

      expect(asked).toEqual([{ dimension: 'toolCalls', need: 1 }]);
      const order = runner.observations;
      const opened = order.findIndex((o) => o.kind === 'decision.opened' && o.decisionKind === 'budget');
      const startedOf = (tool: string) =>
        order.findIndex((o) => o.kind === 'attempt.started' && o.callId === s.tools.callsOf(tool)[0]?.callId);
      expect(opened).toBeGreaterThanOrEqual(0);
      expect(startedOf('a')).toBeLessThan(opened);
      expect(startedOf('b')).toBeGreaterThan(opened);
      expect(s.tools.calls.map((call) => call.toolId)).toEqual(['a', 'b']);
      expect(await statuses(s)).toEqual(['completed']);
    });

    it('V6 price-cost: a price function charges cost on every attempt and a cost limit opens budget (harness.md §11)', async () => {
      // Cost in tokens: input plus ten per output token.
      const price = (_configuration: unknown, usage: Usage) => usage.uncachedInput + 10 * usage.output;
      const asked: { dimension: string; need: number | undefined }[] = [];
      const s = await createSession(fixture, {
        configuration: { limits: { cost: 1_050 } },
        tools: { echo: () => ({ text: 'echo' }) },
        model: [
          { calls: [{ tool: 'echo', arguments: {} }], usage: { inputTokens: 100, outputTokens: 90 } },
          { text: 'final', usage: { inputTokens: 50, outputTokens: 5 } },
        ],
        drive: with_({
          budget: ({ payload }) => {
            asked.push({ dimension: payload.dimension, need: payload.need });
            return { extend: { cost: 100_000 } };
          },
        }),
      });
      const runner = s.runner({ price });
      await s.send('go', { runner });

      const costs = stepEvents(await s.events()).flatMap((event) =>
        event.calls.flatMap((call) => (call.kind === 'model' ? call.attempts.map((a) => a.usage.cost) : [])),
      );
      expect(costs).toEqual([1_000, 100]);
      expect((await s.turns())[0]?.usage.cost).toBe(1_100);
      // The first call fit; the second needed its input estimate in cost on top of the 1,000 spent.
      expect(asked).toHaveLength(1);
      expect(asked[0]?.dimension).toBe('cost');
      expect(asked[0]?.need).toBe(s.model.invocations[1]?.inputTokens);
      const order = runner.observations;
      const opened = order.findIndex((o) => o.kind === 'decision.opened' && o.decisionKind === 'budget');
      // The model call and the tool call started before the budget opened; the second model call after.
      expect(opened).toBeGreaterThanOrEqual(0);
      expect(order.slice(0, opened).filter((o) => o.kind === 'attempt.started')).toHaveLength(2);
      expect(order.slice(opened).some((o) => o.kind === 'attempt.started')).toBe(true);
      expect(await statuses(s)).toEqual(['completed']);
    });

    it("V6 default-budget: defaultDrive's budget without an extension ends the turn failed (harness.md §7, §11)", async () => {
      const s = await createSession(fixture, {
        configuration: { limits: { modelCalls: 1 } },
        tools: { echo: () => ({ text: 'echo' }) },
        model: [{ calls: [{ tool: 'echo', arguments: {} }] }, { text: 'never' }],
      });
      await s.send('go');

      expect(await statuses(s)).toEqual(['failed']);
      const last = lastStep(await s.events());
      expect(last).toMatchObject({ stepKind: 'model', ending: 'abandoned', end: { status: 'failed' } });
      expect(JSON.stringify(last.end?.reason)).toContain('budget');
      expect(s.model.invocations).toHaveLength(1);
      expect(s.tools.callsOf('echo')).toHaveLength(1);
    });

    // V7. Defaults and the schedule (harness.md §6, §7).

    it('V7 max-concurrency: maxConcurrency 1 within a stage of three calls never runs two at once (harness.md §6)', async () => {
      const serial = await concurrency(fixture, 1);
      expect(serial.peak).toBe(1);
      expect(serial.order).toEqual(['a', 'b', 'c']);
      // Without the cap the same stage does run them together.
      expect((await concurrency(fixture, 3)).peak).toBe(3);
    });

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
        const s = await createSession(fixture, setup);
        await s.send('go');
        expectWith(await statuses(s), name).toEqual(['failed']);
        expectWith(JSON.stringify((await endOf(s))?.reason), name).toContain('invalid-answer');
        expectWith(s.tools.calls, name).toHaveLength(0);
      }
    });

    it('V7 default-model-retry: defaultDrive retries a transient model failure three attempts with growing waits, then fails the turn (harness.md §7)', async () => {
      const s = await createSession(fixture, {
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
      const s = await createSession(fixture, {
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

    it('V7 stop-at-permission: stop is accepted at permission and the turn ends stopped (harness.md §7)', async () => {
      const s = await createSession(fixture, {
        configuration: { tools: { enabled: ['t'], preauthorized: [] } },
        tools: { t: () => ({ text: 'ran' }) },
        model: [{ calls: [{ tool: 't', arguments: {} }] }, { text: 'never' }],
        drive: with_({ permission: () => 'stop' }),
      });
      await s.send('go');

      expect(await statuses(s)).toEqual(['stopped']);
      expect(lastStep(await s.events())).toMatchObject({ ending: 'abandoned', end: { status: 'stopped' } });
      expect(s.tools.calls).toHaveLength(0);
      expect(s.model.remaining()).toBe(1);
    });

    for (const order of [
      ['park', 'fail'],
      ['fail', 'park'],
    ] as const) {
      it(`V7 fail-beats-park: a step answered ${order.join(' then ')} ends failed, not parked (harness.md §6)`, async () => {
        const answers = {
          park: { park: { note: 'ask the person' } },
          fail: { fail: 'permission refused' },
        } as const;
        const s = await createSession(fixture, {
          configuration: { tools: { enabled: ['a', 'b'], preauthorized: [] } },
          tools: { a: () => ({ text: 'a' }), b: () => ({ text: 'b' }) },
          model: [
            {
              calls: [
                { tool: 'a', arguments: {} },
                { tool: 'b', arguments: {} },
              ],
            },
          ],
          drive: with_({ permission: ({ payload }) => answers[payload.call.toolId === 'a' ? order[0] : order[1]] }),
        });
        await s.send('go');

        expect(await statuses(s)).toEqual(['failed']);
        const tools = stepEvents(await s.events()).filter((event) => event.stepKind === 'tools');
        expect(tools.map((event) => event.ending)).toEqual(['abandoned']);
        expect(JSON.stringify(tools[0]?.end?.reason)).toContain('permission refused');
        expect((await s.client.execution()).turn).toBeUndefined();
        expect(s.tools.calls).toHaveLength(0);
      });
    }

    // V8. The model result (harness.md §6).

    it('V8 accept-output-limit: accept on an output-limit message with tool calls plans and runs them (harness.md §6, §7)', async () => {
      const seen: string[] = [];
      const s = await createSession(fixture, {
        tools: { echo: () => ({ text: 'echo' }) },
        model: [{ calls: [{ tool: 'echo', arguments: {} }], finish: 'output-limit' }, { text: 'final' }],
        drive: with_({
          'model-result': ({ payload }) => {
            seen.push(`${payload.finish}:${payload.message.parts.filter((p) => p.kind === 'tool-call').length}`);
            return 'accept';
          },
        }),
      });
      await s.send('go');

      expect(seen).toEqual(['output-limit:1']);
      expect(s.tools.callsOf('echo')).toHaveLength(1);
      expect(await resultsOf(s)).toEqual([{ status: 'succeeded', text: 'echo' }]);
      const steps = stepEvents(await s.events());
      const model = steps.find((event) => event.stepKind === 'model');
      expect(model?.next?.kind).toBe('tools');
      expect(steps.filter((event) => event.stepKind === 'tools')).toHaveLength(1);
      expect(await statuses(s)).toEqual(['completed']);
    });

    it('V8 filtered: a filtered finish opens model-result too, and accept completes the turn (harness.md §6)', async () => {
      const seen: string[] = [];
      const s = await createSession(fixture, {
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

    // V9. Interrupts (harness.md §10).

    it('V9 interrupt-idle: an interrupt item naming a turn that is not active is an ordinary item (harness.md §10)', async () => {
      const asked: string[] = [];
      const s = await createSession(fixture, {
        model: [{ text: 'hello' }],
        drive: with_({
          interrupt: ({ payload }) => {
            asked.push(payload.item.key);
            return 'abandon';
          },
        }),
      });
      const { launch } = await s.send('go', { key: 'int', interrupt: { turnId: 'turn-ghost' }, run: false });
      expect(launch).toBeDefined();
      await s.run(launch as string);

      expect(asked).toEqual([]);
      expect((await s.main()).map(textOf)).toEqual(['go', 'hello']);
      expect(await statuses(s)).toEqual(['completed']);
      expect(stepEvents(await s.events()).flatMap((event) => event.consumed)).toEqual(['int']);
    });

    it('V9 interrupt-other-turn: an interrupt item naming another turn does not cut the active one; it begins the next (harness.md §10)', async () => {
      const asked: string[] = [];
      const s = await createSession(fixture, {
        model: [{ text: 'slow', waitMs: 60_000 }, { text: 'second' }],
        drive: with_({
          interrupt: ({ payload }) => {
            asked.push(payload.item.key);
            return 'abandon';
          },
          boundary: () => 'proceed',
        }),
      });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);
      await s.client.submit({ key: 'int', messages: [user('elsewhere')], interrupt: { turnId: 'turn-ghost' } });
      await s.clock.runUntil(running);

      expect(asked).toEqual([]);
      expect(stepEvents(await s.events()).filter((event) => event.ending === 'interrupted')).toHaveLength(0);
      expect(await statuses(s)).toEqual(['completed', 'completed']);
      expect((await s.turns())[0]?.startedBy).toEqual(['int']);
      expect((await s.main()).map(textOf)).toEqual(['go', 'slow', 'elsewhere', 'second']);
    });

    it('V9 interrupt-first: the interrupting item comes first among the steered items in main (harness.md §10)', async () => {
      const later: Later = {};
      const s = await createSession(fixture, {
        tools: {
          a: async (context) => {
            await later.s?.client.submit({ key: 'other', messages: [user('other')] });
            await later.s?.client.submit({
              key: 'int',
              messages: [user('urgent')],
              interrupt: { turnId: context.turnId },
            });
            await nap(later, 1_000, context.signal);
            return { text: 'a done' };
          },
        },
        model: [{ calls: [{ tool: 'a', arguments: {} }] }, { text: 'after' }],
        drive: with_({
          interrupt: () => 'skip',
          // The drive names the ordinary item first; the interrupting one still leads.
          boundary: (_decision, context) =>
            context.inbox.length > 0 ? { steer: [...context.inbox].reverse().map((item) => item.key) } : 'proceed',
        }),
      });
      later.s = s;
      await s.send('go');

      const users = (await s.main()).filter((message) => message.role === 'user').map(textOf);
      expect(users.filter((text) => text !== '')).toEqual(['go', 'urgent', 'other']);
      const settling = stepEvents(await s.events()).find(
        (event) => event.consumed.length > 0 && event.stepKind === 'tools',
      );
      expect(settling?.consumed).toEqual(['int', 'other']);
      expect(await statuses(s)).toEqual(['completed']);
    });

    // V10. Bookkeeping (harness.md §3, §8, §13).

    it("V10 id-prefixes: call and attempt ids are prefixed with their step's id, nested calls included (harness.md §3)", async () => {
      const s = await createSession(fixture, {
        declare: ['outer'],
        tools: {
          outer: async (context) => ({ text: (await context.run(nestedRun('inner'))).status }),
          inner: () => ({ text: 'inner done' }),
        },
        model: [
          { calls: [{ tool: 'outer', arguments: {} }] },
          { fail: { transient: true } },
          { calls: [{ tool: 'inner', arguments: {} }] },
          { text: 'nested done' },
          { text: 'final' },
        ],
      });
      await s.send('go');

      const steps = stepEvents(await s.events()).filter((event) => event.calls.length > 0);
      const calls = steps.flatMap((event) => event.calls);
      expect(steps.length).toBeGreaterThanOrEqual(3);
      expect(calls.filter((call) => call.parentCallId !== undefined).length).toBeGreaterThanOrEqual(3);
      expect(calls.flatMap((call) => call.attempts).length).toBeGreaterThanOrEqual(calls.length + 1);
      for (const event of steps) {
        for (const call of event.calls) {
          expectWith(call.id.startsWith(event.step), `call ${call.id} of ${event.step}`).toBe(true);
          for (const attempt of call.attempts)
            expectWith(attempt.id.startsWith(event.step), `attempt ${attempt.id} of ${event.step}`).toBe(true);
        }
      }
      const attemptIds = calls.flatMap((call) => call.attempts.map((attempt) => attempt.id));
      expect(new Set(idsOf(calls)).size).toBe(calls.length);
      expect(new Set(attemptIds).size).toBe(attemptIds.length);
    });

    it("V10 request-evidence: a model call's request blob holds the adapter revision and the payload digest (harness.md §3)", async () => {
      const s = await createSession(fixture, {
        tools: { echo: () => ({ text: 'echo' }) },
        model: [{ calls: [{ tool: 'echo', arguments: {} }] }, { text: 'final' }],
      });
      await s.send('go');

      const modelCalls = stepEvents(await s.events()).flatMap((event) => event.calls.filter((c) => c.kind === 'model'));
      expect(modelCalls).toHaveLength(2);
      const bodies = (await s.client.blobs(modelCalls.map((call) => call.request))) as Record<string, Json>[];
      expect(bodies.map((body) => body.adapterRevision)).toEqual(['scripted-1', 'scripted-1']);
      expect(bodies.map((body) => body.digest)).toEqual(s.model.invocations.map((invocation) => invocation.digest));
      expect(bodies[0]?.digest).not.toBe(bodies[1]?.digest);
      // What reproduces the call, not its payload.
      expect(bodies[0]?.configuration).toEqual(TURN.model);
      expect(bodies.every((body) => !('payload' in body))).toBe(true);
    });

    it("V10 options-reach-configure: an item's options reach configure, and its configuration is the turn's (harness.md §5)", async () => {
      const seen: (Json | undefined)[][] = [];
      const configure = (items: readonly Item[]) => {
        seen.push(items.map((item) => item.options));
        const options = items.at(-1)?.options as { temperature?: number } | undefined;
        return {
          ...TURN,
          model: {
            ...TURN.model,
            parameters: { ...(options?.temperature !== undefined ? { temperature: options.temperature } : {}) },
          },
        };
      };
      const s = await createSession(fixture, { model: [{ text: 'hello' }] });
      const { launch } = await s.send('one', { key: 'one', options: { n: 1 }, run: false });
      await s.client.submit({ key: 'two', messages: [user('two')], options: { n: 2, temperature: 0.25 } });
      await s.run(launch as string, s.runner({ configure }));

      expect(seen.at(-1)).toEqual([{ n: 1 }, { n: 2, temperature: 0.25 }]);
      const begin = stepEvents(await s.events()).find((event) => event.stepKind === 'begin');
      expect(begin?.begun?.items).toEqual(['one', 'two']);
      expect(begin?.begun?.configuration.model.parameters.temperature).toBe(0.25);
      expect(s.model.projections[0]?.configuration.parameters.temperature).toBe(0.25);
      expect(await statuses(s)).toEqual(['completed']);
    });

    it("V10 last-usage-frame: a handler's last usage frame is what the attempt records (harness.md §8)", async () => {
      const s = await createSession(fixture, {
        tools: {
          metered: async function* (): AsyncGenerator<ToolFrame> {
            yield { kind: 'usage', usage: { uncachedInput: 5 } };
            yield { kind: 'progress', value: 'working' };
            yield { kind: 'usage', usage: { uncachedInput: 12, output: 3 } };
            yield { kind: 'result', status: 'succeeded', parts: [{ kind: 'text', text: 'done' }] };
          },
        },
        model: [{ calls: [{ tool: 'metered', arguments: {} }] }, { text: 'final' }],
      });
      await s.send('go');

      const call = stepEvents(await s.events())
        .flatMap((event) => event.calls)
        .find((c) => c.kind === 'tool');
      expect(call?.attempts).toHaveLength(1);
      expect(call?.attempts[0]?.usage).toMatchObject({ uncachedInput: 12, output: 3, toolCalls: 1 });
    });

    it("V10 journal-edit: journal.edit on a pointer is staged and moves it with the step's commit (harness.md §8)", async () => {
      const later: Later = {};
      const [first, second, third] = [user('first'), user('second'), user('third')];
      const seen: { before?: string[]; during?: string[]; stored?: string | undefined; head?: string | undefined } = {};
      const s = await createSession(fixture, {
        tools: {
          seed: (context) => {
            context.journal.append('log', [first, second, third]);
            return { text: 'seeded' };
          },
          squash: async (context) => {
            const before = await context.journal.read('log');
            seen.before = before.messages.map(textOf);
            seen.head = before.head;
            context.journal.edit('log', { entries: [await hashMessage(first), user('summary')] });
            seen.during = (await context.journal.read('log')).messages.map(textOf);
            seen.stored = (await later.s?.client.pointers(['log']))?.[0]?.head;
            return { text: 'squashed' };
          },
        },
        model: [
          { calls: [{ tool: 'seed', arguments: {} }] },
          { calls: [{ tool: 'squash', arguments: {} }] },
          { text: 'final' },
        ],
      });
      later.s = s;
      await s.send('go');

      expect(seen.before).toEqual(['first', 'second', 'third']);
      expect(seen.during).toEqual(['first', 'summary']);
      // Inside the step the pointer had not moved.
      expect(seen.head).toBeDefined();
      expect(seen.stored).toBe(seen.head);
      const [pointer] = await s.client.pointers(['log']);
      expect(pointer?.head).not.toBe(seen.head);
      expect((await messagesOf(s, 'log')).map(textOf)).toEqual(['first', 'summary']);
      const step = stepEvents(await s.events()).find((event) => event.moves.some((move) => move.from === seen.head));
      expect(step?.moves.find((move) => move.pointer === 'log')).toMatchObject({ from: seen.head, to: pointer?.head });
      expect(step?.compaction?.map((c) => c.pointer)).toEqual(['log']);
      expect(step?.moves.map((move) => move.pointer).sort()).toEqual(['log', MAIN]);
    });

    it('V10 observations: model.delta, attempt.ended and decision.answered arrive in order and survive JSON (harness.md §13)', async () => {
      const s = await createSession(fixture, {
        model: [{ fail: { transient: true } }, { reasoning: 'think', deltas: ['he', 'llo'], text: 'hello' }],
      });
      const runner = s.runner();
      await s.send('go', { runner });

      const live = runner.observations.filter((o): o is Exclude<Observation, { kind: 'event' }> => o.kind !== 'event');
      for (const observation of live) expect(JSON.parse(JSON.stringify(observation))).toEqual(observation);

      const deltas = live.flatMap((o) => (o.kind === 'model.delta' ? [o.delta] : []));
      expect(deltas).toEqual([{ reasoning: 'think' }, { text: 'he' }, { text: 'llo' }]);
      const ended = live.flatMap((o) => (o.kind === 'attempt.ended' ? [o] : []));
      expect(ended.map((o) => o.outcome)).toEqual(['failed', 'succeeded']);
      expect(ended[0]?.callId).toBe(ended[1]?.callId);
      expect(ended[0]?.attemptId).not.toBe(ended[1]?.attemptId);
      for (const o of ended) {
        expect(o.path).toEqual([]);
        expect(o.attemptId.startsWith(o.callId)).toBe(true);
      }
      const answered = live.flatMap((o) => (o.kind === 'decision.answered' ? [`${o.decisionKind}:${o.answer}`] : []));
      expect(answered).toEqual(['next-turn:begin', 'failure:retry', 'boundary:proceed', 'next-turn:release']);
      for (const o of live)
        if (o.kind === 'decision.answered') {
          expect(o.path).toEqual([]);
          expect(o.ids).toBeDefined();
        }
      // Order within the second attempt: started, deltas, ended; then the boundary opens and is answered.
      const at = (match: (o: Exclude<Observation, { kind: 'event' }>) => boolean) => live.findIndex(match);
      const started = at((o) => o.kind === 'attempt.started' && o.attemptId === ended[1]?.attemptId);
      const delta = at((o) => o.kind === 'model.delta');
      const done = at((o) => o.kind === 'attempt.ended' && o.attemptId === ended[1]?.attemptId);
      const opened = at((o) => o.kind === 'decision.opened' && o.decisionKind === 'boundary');
      const closed = at((o) => o.kind === 'decision.answered' && o.decisionKind === 'boundary');
      expect(started).toBeLessThan(delta);
      expect(delta).toBeLessThan(done);
      expect(done).toBeLessThan(opened);
      expect(opened).toBeLessThan(closed);
    });
  });
}
