/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { digest, hashMessage, type Json, type Message } from '@molstar/ai-agent';
import { type CallRecord, type CommitBatch, MAIN, SCHEMA_VERSION, type Usage } from '@molstar/ai-harness';
import { type ConformanceFixture, createSession, type Session, type SessionSetup } from '../session.js';
import { settle } from '../testing/clock.js';
import { assistant, system, textOf, user } from '../testing/messages.js';
import { appendMove, expectWith, stepEvents, TURN, with_ } from './helpers.js';

/** Waits (in real turns of the event loop) until `ok` holds. */
async function until(ok: () => boolean | Promise<boolean>, rounds = 200): Promise<void> {
  for (let i = 0; i < rounds; i += 1) {
    if (await ok()) return;
    await settle();
  }
  throw new Error('until: the condition never held');
}

const NO_USAGE: Usage = {
  modelCalls: 0,
  toolCalls: 0,
  uncachedInput: 0,
  cacheRead: 0,
  cacheWrite: 0,
  output: 0,
  reasoning: 0,
};

/** Protocol level: claims and begins a turn over item "item-1". The model step "step-2" is pending. */
async function beginTurn(s: Session) {
  const { launch } = await s.client.submit({ key: 'item-1', messages: [user('hello')] });
  const { lease, execution } = await s.port.claim({
    launch: launch as string,
    runnerId: 'r1',
    schemaVersion: SCHEMA_VERSION,
  });
  const begun = await s.port.begin(lease, {
    key: 'begin-1',
    turnId: 'turn-x',
    step: 'step-1',
    items: ['item-1'],
    configuration: TURN,
    messages: [user('hello')],
    move: await appendMove(s.port, MAIN, execution.main, [user('hello')]),
    calls: [],
    next: { id: 'step-2', kind: 'model' },
  });
  return begun.lease;
}

/** Protocol level: as `beginTurn`, then the model step settles with one call; the tools step "step-3" is pending. */
async function toolsPending(s: Session) {
  const lease = await beginTurn(s);
  const message: Message = {
    role: 'assistant',
    parts: [{ kind: 'tool-call', requestId: 'req-1', tool: { id: 't' }, arguments: {} }],
  };
  const body = message as unknown as Json;
  const evidence: Json = { request: 'evidence-1' };
  const messageRef = await digest(body);
  const requestRef = await digest(evidence);
  await s.port.putBlobs(lease, [
    { ref: messageRef, body },
    { ref: requestRef, body: evidence },
  ]);
  const now = s.clock.now();
  const call: CallRecord = {
    id: 'call-1',
    key: 'call-1',
    kind: 'model',
    request: requestRef,
    attempts: [
      {
        id: 'attempt-1',
        startedAt: now,
        endedAt: now,
        outcome: 'succeeded',
        usage: { ...NO_USAGE, modelCalls: 1 },
      },
    ],
  };
  const committed = await s.port.commit(lease, {
    key: 'commit-model',
    turnId: 'turn-x',
    step: 'step-2',
    ending: 'settled',
    calls: [call],
    messages: [],
    moves: [],
    consumed: [],
    next: {
      id: 'step-3',
      kind: 'tools',
      message: messageRef,
      plan: [{ callId: 'call-2', key: 'call-2', requestId: 'req-1', toolId: 't', arguments: {} }],
    },
  });
  return committed.lease;
}

/**
 * A top-level tool `outer` runs a nested run whose first round calls `inner` (preauthorized) and
 * whose second calls `guarded`, whose `permission` the app's drive answers `park`.
 */
async function nestedPark(fixture: ConformanceFixture) {
  const asked: { path: readonly string[]; toolId: string }[] = [];
  let session: Session | undefined;
  const s = await createSession(fixture, {
    tools: {
      outer: async (context) => {
        const result = await context.run({
          input: [user('investigate')],
          tools: {
            enabled: ['inner', 'guarded'],
            preauthorized: ['inner'],
            definitions: [session?.tools.definition('inner'), session?.tools.definition('guarded')].flatMap((d) =>
              d === undefined ? [] : [d],
            ),
          },
          model: TURN.model,
          limits: {},
        });
        return { text: result.status };
      },
      inner: () => ({ text: 'inner done' }),
      guarded: () => ({ text: 'guarded done' }),
    },
    declare: ['outer'],
    configuration: { tools: { enabled: ['outer'], preauthorized: ['outer'] } },
    model: [
      { calls: [{ tool: 'outer', arguments: {} }] },
      { calls: [{ tool: 'inner', arguments: {} }], usage: { inputTokens: 700, outputTokens: 30 } },
      { calls: [{ tool: 'guarded', arguments: {} }], usage: { inputTokens: 700, outputTokens: 30 } },
    ],
    drive: with_({
      permission: ({ path, payload }) => {
        asked.push({ path, toolId: payload.call.toolId });
        return payload.call.toolId === 'guarded' ? { park: { note: 'ask the person' } } : 'allow';
      },
    }),
  });
  session = s;
  return { s, asked };
}

/** Steps (conformance.md): E1–E6. */
export function stepCases(fixture: ConformanceFixture): void {
  describe('Steps', () => {
    // E1. A step ends exactly one way, its commit matches that ending, and nothing is decided after the commit.

    it('E1 parked-with-move: a parked batch with a move fails invalid', async () => {
      const s = await createSession(fixture);
      const lease = await toolsPending(s);
      const main = (await s.client.execution()).main;
      const move = await appendMove(s.port, MAIN, main, [assistant('x')]);
      const parked = (over: Partial<CommitBatch>): CommitBatch => ({
        key: 'park',
        turnId: 'turn-x',
        step: 'step-3',
        ending: 'parked',
        calls: [],
        messages: [],
        moves: [],
        consumed: [],
        park: { path: [] },
        ...over,
      });
      // harness.md §6: a park moves nothing.
      await expect(
        s.port.commit(lease, parked({ key: 'park-bad', messages: [assistant('x')], moves: [move] })),
      ).rejects.toMatchObject({ code: 'invalid' });
      const untouched = await s.client.execution();
      expect(untouched.main).toBe(main);
      expect(untouched.turn?.park).toBeUndefined();
      // The same commit without the move parks.
      await s.port.commit(lease, parked({ key: 'park-ok' }));
      expect((await s.client.execution()).turn?.park).toBeDefined();
    });

    it('E1 abandoned-with-consumed: an abandoned batch consuming an item fails invalid', async () => {
      const s = await createSession(fixture);
      const lease = await beginTurn(s);
      await s.client.submit({ key: 'extra', messages: [user('extra')] });
      const abandoned = (over: Partial<CommitBatch>): CommitBatch => ({
        key: 'abandon',
        turnId: 'turn-x',
        step: 'step-2',
        ending: 'abandoned',
        calls: [],
        messages: [],
        moves: [],
        consumed: [],
        end: { status: 'failed' },
        ...over,
      });
      // harness.md §6: an abandoned step consumes nothing.
      await expect(s.port.commit(lease, abandoned({ key: 'abandon-bad', consumed: ['extra'] }))).rejects.toMatchObject({
        code: 'invalid',
      });
      expect((await s.client.execution()).turn).toBeDefined();
      await s.port.commit(lease, abandoned({ key: 'abandon-ok' }));
      const after = await s.client.execution();
      expect(after.turn).toBeUndefined();
      expect(after.inbox.map((item) => item.key)).toEqual(['extra']);
    });

    it('E1 interrupted-shape: an interrupted batch appends only steered items and names a model step', async () => {
      const s = await createSession(fixture);
      const lease = await beginTurn(s);
      await s.client.submit({ key: 'int-1', messages: [user('redirect')], interrupt: { turnId: 'turn-x' } });
      const main = (await s.client.execution()).main;
      const interrupted = (over: Partial<CommitBatch>): CommitBatch => ({
        key: 'interrupt',
        turnId: 'turn-x',
        step: 'step-2',
        ending: 'interrupted',
        calls: [],
        messages: [user('redirect')],
        moves: [],
        consumed: ['int-1'],
        next: { id: 'step-3', kind: 'model' },
        ...over,
      });
      const steered = await appendMove(s.port, MAIN, main, [user('redirect')]);
      const extra = await appendMove(s.port, MAIN, main, [user('redirect'), assistant('extra')]);
      // More than the steered items.
      await expect(
        s.port.commit(
          lease,
          interrupted({ key: 'int-bad-1', messages: [user('redirect'), assistant('extra')], moves: [extra] }),
        ),
      ).rejects.toMatchObject({ code: 'invalid' });
      // Names no model step.
      await expect(
        s.port.commit(lease, interrupted({ key: 'int-bad-2', moves: [steered], end: { status: 'completed' } })),
      ).rejects.toMatchObject({ code: 'invalid' });
      const { next: _next, ...withoutNext } = interrupted({ key: 'int-bad-3', moves: [steered] });
      await expect(s.port.commit(lease, withoutNext)).rejects.toMatchObject({ code: 'invalid' });
      expect((await s.client.execution()).main).toBe(main);
      // The legal shape is accepted.
      await s.port.commit(lease, interrupted({ key: 'int-ok', moves: [steered] }));
      const last = stepEvents(await s.events()).at(-1);
      expect(last?.ending).toBe('interrupted');
      expect(last?.consumed).toEqual(['int-1']);
      expect(last?.next?.kind).toBe('model');
      expect((await s.client.execution()).inbox).toEqual([]);
    });

    it('E1 crash-before-decision: a crash while model-result is open re-runs the step; one settled call (defect 4)', async () => {
      let open = false;
      const s = await createSession(fixture, { model: [{ text: 'cut', finish: 'output-limit' }, { text: 'full' }] });
      const { launch } = await s.send('go', { run: false });
      const doomed = s.runner({
        drive: with_({
          'model-result': () => {
            open = true;
            return new Promise<never>(() => undefined);
          },
        }),
      });
      doomed.run(launch as string).catch(() => undefined);
      await until(() => open);
      // harness.md §6: a crash records nothing, so the open decision has left no model step event.
      expect(stepEvents(await s.events()).filter((e) => e.stepKind === 'model')).toHaveLength(0);
      doomed.abort('crash');
      await settle();
      await s.clock.advance(3_600_000);
      const { launch: relaunch } = await s.client.check();
      expect(relaunch).toBeDefined();
      await s.run(relaunch as string, s.runner({ drive: with_({ 'model-result': () => 'accept' }) }));

      const modelSteps = stepEvents(await s.events()).filter((e) => e.stepKind === 'model');
      expect(modelSteps).toHaveLength(1);
      expect(modelSteps[0]?.ending).toBe('settled');
      const calls = modelSteps[0]?.calls ?? [];
      expect(calls).toHaveLength(1);
      expect(calls[0]?.attempts.map((a) => a.outcome)).toEqual(['succeeded']);
      const [turn] = await s.turns();
      expect(turn?.status).toBe('completed');
      expect((await s.main()).map((m) => textOf(m))).toContain('full');
    });

    it('E1 invalid-answer: an answer not listed for a kind fails the run with invalid-answer, at every kind (defect 14)', async () => {
      const tool = { t: () => ({ text: 'done' }) };
      const call = { calls: [{ tool: 't', arguments: {} }] };
      const cases: [string, SessionSetup][] = [
        [
          'permission',
          {
            tools: tool,
            configuration: { tools: { enabled: ['t'], preauthorized: [] } },
            model: [call],
            drive: with_({ permission: () => ({ steer: ['item-x'] }) as never }),
          },
        ],
        ['schedule', { tools: tool, model: [call], drive: with_({ schedule: () => 'allow' as never }) }],
        [
          'budget',
          {
            tools: tool,
            configuration: { limits: { modelCalls: 1 } },
            model: [call, { text: 'x' }],
            drive: with_({ budget: () => 'proceed' as never }),
          },
        ],
        ['failure', { model: [{ fail: {} }], drive: with_({ failure: () => 'allow' as never }) }],
        [
          'model-result',
          {
            model: [{ text: 'cut', finish: 'output-limit' }],
            drive: with_({ 'model-result': () => 'proceed' as never }),
          },
        ],
        [
          'repair',
          {
            tools: tool,
            model: [{ calls: [{ tool: 'ghost', arguments: {} }] }],
            drive: with_({ repair: () => 'allow' as never }),
          },
        ],
        ['boundary', { model: [{ text: 'hi' }], drive: with_({ boundary: () => ({ retry: {} }) as never }) }],
      ];
      for (const [kind, setup] of cases) {
        const s = await createSession(fixture, setup);
        await s.send('go');
        const [turn] = await s.turns();
        expectWith(turn?.status, kind).toBe('failed');
        const end = stepEvents(await s.events()).find((e) => e.end !== undefined)?.end;
        expectWith(end?.status, kind).toBe('failed');
        // harness.md §7: any other answer is `fail` with invalid-answer.
        expectWith(JSON.stringify(end?.reason), kind).toContain('invalid-answer');
        if (kind === 'permission') expectWith(s.tools.calls, kind).toHaveLength(0);
      }
      // interrupt: opened by an interrupt item during a model step.
      const s = await createSession(fixture, {
        model: [{ text: 'slow', waitMs: 60_000 }],
        drive: with_({ interrupt: () => 'proceed' as never }),
      });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);
      const turnId = (await s.client.execution()).turn?.id as string;
      await s.client.submit({ key: 'int-1', messages: [user('redirect')], interrupt: { turnId } });
      await s.clock.runUntil(running);
      const end = stepEvents(await s.events()).find((e) => e.end !== undefined)?.end;
      expectWith(end?.status, 'interrupt').toBe('failed');
      expectWith(JSON.stringify(end?.reason), 'interrupt').toContain('invalid-answer');
    });

    // E2. A round moves main once.

    it('E2 round-order: calls finishing in reverse order are appended in call order', async () => {
      let session: Session | undefined;
      const slow = (ms: number, name: string) => async () => {
        await session?.clock.sleep(ms);
        // The fastest call submits the item that will be steered in at the round's boundary.
        if (name === 'c') await session?.client.submit({ key: 'steer-1', messages: [user('steered')] });
        return { text: `after ${ms}`, contributions: [system(`note ${name}`)] };
      };
      const s = await createSession(fixture, {
        tools: { a: slow(300, 'a'), b: slow(200, 'b'), c: slow(100, 'c') },
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
          schedule: ({ payload }) => ({ stages: [payload.calls.map((call) => call.callId)], maxConcurrency: 3 }),
          boundary: ({ payload }, ctx) =>
            payload.completing || ctx.inbox.length === 0 ? 'proceed' : { steer: ctx.inbox.map((item) => item.key) },
        }),
      });
      session = s;
      await s.send('go');
      const main = await s.main();
      const callIndex = main.findIndex((m) => m.role === 'assistant' && m.parts.some((p) => p.kind === 'tool-call'));
      const calls = main[callIndex] as Message;
      const results = main[callIndex + 1] as Message;
      const callIds = calls.parts.flatMap((p) => (p.kind === 'tool-call' ? [p.requestId] : []));
      const resultIds = results.parts.flatMap((p) => (p.kind === 'tool-result' ? [p.requestId] : []));
      expect(results.role).toBe('user');
      expect(resultIds).toEqual(callIds);
      // Then each call's contributions in call order, then the steered item, then the next model step.
      const after = main.slice(callIndex + 2).map((m) => textOf(m));
      expect(after).toEqual(['note a', 'note b', 'note c', 'steered', 'final']);
    });

    it('E2 no-calls-round: a model step without calls appends its message alone and ends the turn', async () => {
      const s = await createSession(fixture, { model: [{ text: 'hello' }] });
      await s.send('go');
      expect((await s.main()).map((m) => m.role)).toEqual(['user', 'assistant']);
      const last = stepEvents(await s.events()).at(-1);
      expect(last?.stepKind).toBe('model');
      expect(last?.moves).toHaveLength(1);
      expect(last?.end?.status).toBe('completed');
      expect(last?.next).toBeUndefined();
    });

    it('E2 skill-twice: two tools contributing the same skill text in one round leave one entry', async () => {
      const skill = system('skill: how to look closely');
      const s = await createSession(fixture, {
        tools: { a: () => ({ text: 'a', contributions: [skill] }), b: () => ({ text: 'b', contributions: [skill] }) },
        model: [
          {
            calls: [
              { tool: 'a', arguments: {} },
              { tool: 'b', arguments: {} },
            ],
          },
          { text: 'final' },
        ],
      });
      await s.send('go');
      const main = await s.main();
      expect(main.filter((m) => textOf(m) === textOf(skill))).toHaveLength(1);
    });

    it('E2 same-input-twice: two identical user messages in a row are both appended', async () => {
      const s = await createSession(fixture, { model: [{ text: 'ok' }] });
      const { launch } = await s.send('same', { run: false });
      await s.send('same', { key: 'item-2', run: false });
      await s.run(launch as string);
      const main = await s.main();
      expect(main.map((m) => (m.role === 'user' ? textOf(m) : m.role))).toEqual(['same', 'same', 'assistant']);
      // Two positions of the head, one stored entry.
      const hash = await hashMessage(user('same'));
      const execution = await s.client.execution();
      const [head] = await s.client.heads([execution.main]);
      expect(head?.entries.filter((entry) => entry === hash)).toHaveLength(2);
      expect(await s.client.entries([hash])).toEqual([user('same')]);
    });

    it("E2 after-edit: a compaction edit removes a skill's entry; a later identical contribution enters again", async () => {
      const skill = system('skill: how to look closely');
      const skillHash = await hashMessage(skill);
      let compacted = false;
      const s = await createSession(fixture, {
        tools: { look: () => ({ text: 'looked', contributions: [skill] }) },
        model: [
          { calls: [{ tool: 'look', arguments: {} }] },
          { calls: [{ tool: 'look', arguments: {} }] },
          { text: 'final' },
        ],
        drive: with_({
          boundary: ({ payload }, ctx) => {
            if (payload.completing || compacted) return 'proceed';
            compacted = true;
            return { compact: { entries: ctx.head.entries.filter((entry) => entry !== skillHash) } };
          },
        }),
      });
      await s.send('go');
      const main = await s.main();
      const at = main.map((m, i) => (textOf(m) === textOf(skill) ? i : -1)).filter((i) => i >= 0);
      expect(at).toHaveLength(1);
      // It is the second round's contribution: it follows the last tool-result message.
      const lastResults = main.findLastIndex((m) => m.parts.some((p) => p.kind === 'tool-result'));
      expect(at[0]).toBeGreaterThan(lastResults);
    });

    // E3. A step that does not settle records its ended attempts and usage, and no results.

    it("E3 park-usage: a nested run's ended attempts and usage are in the parked step event", async () => {
      const { s } = await nestedPark(fixture);
      await s.send('go');
      const parked = stepEvents(await s.events()).filter((e) => e.ending === 'parked');
      expect(parked).toHaveLength(1);
      const [event] = parked;
      expect(event?.moves).toEqual([]);
      const nested = (event?.calls ?? []).filter((c) => c.parentCallId !== undefined);
      // The nested tool attempt that ended before the park ...
      const inner = nested.filter((c) => c.kind === 'tool');
      expect(inner).toHaveLength(1);
      expect(inner[0]?.attempts.map((a) => a.outcome)).toEqual(['succeeded']);
      // ... and both nested model attempts with their usage.
      const models = nested.filter((c) => c.kind === 'model');
      expect(models.map((c) => c.attempts[0]?.usage.uncachedInput)).toEqual([700, 700]);
      expect(models.map((c) => c.attempts[0]?.usage.output)).toEqual([30, 30]);
      // The turn's running totals include them (harness.md §3, §4).
      expect(((await s.client.execution()).turn?.usage.uncachedInput ?? 0) >= 1_400).toBe(true);
    });

    it('E3 crash-nothing: a crash mid-step writes no event', async () => {
      const s = await createSession(fixture, { model: [{ text: 'slow', waitMs: 60_000 }] });
      const { launch } = await s.send('go', { run: false });
      const doomed = s.runner();
      doomed.run(launch as string).catch(() => undefined);
      await until(() => s.model.invocations.length === 1);
      await settle();
      const before = await s.events();
      const seq = (await s.client.execution()).seq;
      expect(stepEvents(before).filter((e) => e.stepKind === 'model')).toHaveLength(0);
      doomed.abort('crash');
      await settle();
      await s.clock.advance(5_000);
      expect(await s.events()).toEqual(before);
      expect((await s.client.execution()).seq).toBe(seq);
    });

    // E4. A retry is a new attempt of the same call, with the same request; `replace` is a new call.

    it('E4 retry-same-digest: a failed model attempt and a retry make one call with two attempts of equal digests', async () => {
      const s = await createSession(fixture, { model: [{ fail: { transient: true } }, { text: 'ok' }] });
      await s.send('go');
      const calls = stepEvents(await s.events()).flatMap((e) => e.calls.filter((c) => c.kind === 'model'));
      expect(calls).toHaveLength(1);
      expect(calls[0]?.attempts.map((a) => a.outcome)).toEqual(['failed', 'succeeded']);
      expect(s.model.invocations).toHaveLength(2);
      const [first, second] = s.model.invocations;
      expect(second?.digest).toBe(first?.digest);
      expect(second?.adapterRevision).toBe(first?.adapterRevision);
    });

    it('E4 retry-changed-request: an adapter that projects a different request on retry fails the attempt integrity', async () => {
      let session: Session | undefined;
      const failures: string[] = [];
      const s = await createSession(fixture, {
        model: [{ fail: { transient: true } }, { text: 'ok' }],
        // The payload changes once an attempt has been invoked: the retry's projection differs.
        modelOptions: {
          payload: (input) => ({
            messages: input.messages as unknown as Json,
            invoked: session?.model.invocations.length ?? 0,
          }),
        },
        drive: with_({
          failure: ({ payload }) => {
            failures.push(payload.failure.code);
            return failures.length === 1 ? { retry: {} } : { fail: 'enough' };
          },
        }),
      });
      session = s;
      await s.send('go');
      // harness.md §6: the same digest and adapter revision, or the attempt fails integrity.
      expect(failures.at(-1)).toBe('integrity');
      const call = stepEvents(await s.events())
        .flatMap((e) => e.calls)
        .find((c) => c.kind === 'model');
      const last = call?.attempts.at(-1);
      expect(last?.outcome).toBe('failed');
      expect(last?.failure?.code).toBe('integrity');
    });

    it('E4 retry-backoff: a retry wait longer than signalEveryMs writes no event and keeps the lease live', async () => {
      let retrying = false;
      const s = await createSession(fixture, {
        host: { leaseMs: 3_000 },
        runnerLimits: { signalEveryMs: 1_000 },
        model: [{ fail: { transient: true } }, { text: 'ok' }],
        drive: with_({
          failure: () => {
            retrying = true;
            return { retry: { afterMs: 10_000 } };
          },
        }),
      });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => retrying);
      await settle();
      const before = await s.events();
      const seq = (await s.client.execution()).seq;
      // Twice the lease passes inside the wait: only renewals keep the runner's hold.
      await s.clock.advance(6_000);
      expect(await s.events()).toEqual(before);
      const after = await s.client.execution();
      expect(after.seq).toBe(seq);
      expect(after.runner.phase).toBe('live');
      expect(after.runner.phase === 'live' && after.runner.leaseUntil > s.clock.now()).toBe(true);
      await s.clock.runUntil(running);
      const [turn] = await s.turns();
      expect(turn?.status).toBe('completed');
      expect(s.model.invocations).toHaveLength(2);
    });

    it('E4 replace-new-call: replace after output-limit makes a second call in the same step', async () => {
      const bigger = { ...TURN.model, parameters: { maxOutputTokens: 8_000 } };
      const s = await createSession(fixture, {
        model: [{ text: 'cut', finish: 'output-limit' }, { text: 'complete' }],
        drive: with_({ 'model-result': () => ({ replace: bigger }) }),
      });
      await s.send('go');
      const modelSteps = stepEvents(await s.events()).filter((e) => e.stepKind === 'model');
      expect(modelSteps).toHaveLength(1);
      const calls = modelSteps[0]?.calls ?? [];
      expect(calls).toHaveLength(2);
      expect(new Set(calls.map((c) => c.id)).size).toBe(2);
      expect(calls.map((c) => c.attempts.length)).toEqual([1, 1]);
      // The second call ran with the replacement configuration (for this step only, harness.md §7).
      expect(s.model.projections.at(-1)?.configuration).toEqual(bigger);
      expect(s.model.projections[0]?.configuration).toEqual(TURN.model);
    });

    // E5. Steered input lands only at a round boundary.

    it('E5 no-steer-after-call-step: a model step with calls opens no boundary; the item is offered after the tools step', async () => {
      const timeline: string[] = [];
      const s = await createSession(fixture, {
        tools: {
          work: () => {
            timeline.push('tool');
            return { text: 'worked' };
          },
        },
        model: [{ calls: [{ tool: 'work', arguments: {} }], waitMs: 5_000 }, { text: 'final' }],
        drive: with_({
          boundary: ({ payload }, ctx) => {
            timeline.push(`boundary:${payload.completing ? 'completing' : 'round'}:${ctx.inbox.map((i) => i.key)}`);
            return ctx.inbox.length > 0 && !payload.completing
              ? { steer: ctx.inbox.map((item) => item.key) }
              : 'proceed';
          },
        }),
      });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);
      await s.client.submit({ key: 'late', messages: [user('late')] });
      await s.clock.runUntil(running);
      // No boundary between the model step and the tools step; the first one comes after the tool ran.
      expect(timeline).toEqual(['tool', 'boundary:round:late', 'boundary:completing:']);
      const main = (await s.main()).map((m) => textOf(m));
      expect(main.indexOf('late')).toBeGreaterThan(main.indexOf('worked'));
    });

    it('E5 completing-steer: an item waiting at a completing boundary is steered and the turn continues', async () => {
      const timeline: string[] = [];
      const s = await createSession(fixture, {
        model: [{ text: 'first', waitMs: 5_000 }, { text: 'second' }],
        drive: with_({
          boundary: ({ payload }, ctx) => {
            timeline.push(`boundary:${payload.completing ? 'completing' : 'round'}:${ctx.inbox.map((i) => i.key)}`);
            return ctx.inbox.length > 0 ? { steer: ctx.inbox.map((item) => item.key) } : 'proceed';
          },
        }),
      });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);
      await s.client.submit({ key: 'late', messages: [user('late')] });
      await s.clock.runUntil(running);
      expect(timeline).toEqual(['boundary:completing:late', 'boundary:completing:']);
      expect(s.model.invocations).toHaveLength(2);
      expect((await s.main()).slice(-3).map((m) => `${m.role}:${textOf(m)}`)).toEqual([
        'assistant:first',
        'user:late',
        'assistant:second',
      ]);
      const [turn] = await s.turns();
      expect(turn?.status).toBe('completed');
    });

    it('E5 interrupt-once: an interrupt item cuts one step, is consumed by the commit that follows, and is not seen again', async () => {
      const seen: string[] = [];
      const s = await createSession(fixture, {
        model: [{ text: 'slow', waitMs: 60_000 }, { text: 'after' }],
        drive: with_({
          interrupt: ({ payload }) => {
            seen.push(payload.item.key);
            return 'abandon';
          },
        }),
      });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 1);
      const turnId = (await s.client.execution()).turn?.id as string;
      await s.client.submit({ key: 'int-1', messages: [user('redirect')], interrupt: { turnId } });
      await s.clock.runUntil(running);
      expect(seen).toEqual(['int-1']);
      const steps = stepEvents(await s.events());
      const cut = steps.filter((e) => e.ending === 'interrupted');
      expect(cut).toHaveLength(1);
      expect(cut[0]?.consumed).toEqual(['int-1']);
      expect(cut[0]?.next?.kind).toBe('model');
      expect(steps.filter((e) => e.consumed.includes('int-1'))).toHaveLength(1);
      expect((await s.client.execution()).inbox).toEqual([]);
      const main = (await s.main()).map((m) => `${m.role}:${textOf(m)}`);
      expect(main).toContain('user:redirect');
      expect(main).not.toContain('assistant:slow');
      expect(main.at(-1)).toBe('assistant:after');
    });

    // E6. A park at any depth parks the top-level step; `resume` re-runs it with the drive asked again.

    it('E6 park-resume: permission answered park releases; resume re-runs the step, asks again, settles once', async () => {
      let asked = 0;
      const s = await createSession(fixture, {
        tools: { t: () => ({ text: 'done' }) },
        configuration: { tools: { enabled: ['t'], preauthorized: [] } },
        model: [{ calls: [{ tool: 't', arguments: {} }] }, { text: 'final' }],
        drive: with_({ permission: () => (++asked === 1 ? { park: { note: 'ask the person' } } : 'allow') }),
      });
      const first = await s.send('go');
      expect(first.exit?.reason).toBe('released');
      const parked = await s.client.execution();
      expect(parked.turn?.park).toBeDefined();
      expect(s.tools.calls).toHaveLength(0);

      const { launch } = await s.client.resume({ key: 'resume-1', turnId: parked.turn?.id as string });
      expect(launch).toBeDefined();
      await s.run(launch as string);
      expect(asked).toBe(2);
      expect(s.tools.callsOf('t')).toHaveLength(1);
      const [turn] = await s.turns();
      expect(turn?.status).toBe('completed');
      const toolsStep = stepEvents(await s.events()).find((e) => e.stepKind === 'tools');
      const endings = (await s.client.steps([toolsStep?.step as string])).map((e) => e.ending);
      expect(endings).toEqual(['parked', 'settled']);
    });

    it("E6 nested-park: a nested run's permission answered park parks the top-level step", async () => {
      const { s, asked } = await nestedPark(fixture);
      const { exit } = await s.send('go');
      expect(exit?.reason).toBe('released');
      const guarded = asked.find((a) => a.toolId === 'guarded');
      expect(guarded?.path).toHaveLength(1);
      const toolsStep = stepEvents(await s.events()).find((e) => e.stepKind === 'tools');
      expect(toolsStep?.ending).toBe('parked');
      const execution = await s.client.execution();
      expect(execution.turn?.pending.kind).toBe('tools');
      expect(execution.turn?.park?.path).toEqual(guarded?.path);
      // The top-level step is parked, not the nested run: nothing ran past the park.
      expect(s.tools.callsOf('guarded')).toHaveLength(0);
    });

    it('E6 items-behind-park: items submitted while parked do not launch a runner', async () => {
      const s = await createSession(fixture, {
        tools: { t: () => ({ text: 'done' }) },
        configuration: { tools: { enabled: ['t'], preauthorized: [] } },
        model: [{ calls: [{ tool: 't', arguments: {} }] }],
        drive: with_({ permission: () => ({ park: {} }) }),
      });
      await s.send('go');
      expect((await s.client.execution()).turn?.park).toBeDefined();
      const submitted = await s.client.submit({ key: 'behind', messages: [user('behind the park')] });
      expect(submitted.launch).toBeUndefined();
      expect((await s.client.check()).launch).toBeUndefined();
      const after = await s.client.execution();
      expect(after.runner.phase).toBe('none');
      expect(after.inbox.map((item) => item.key)).toEqual(['behind']);
    });
  });
}
