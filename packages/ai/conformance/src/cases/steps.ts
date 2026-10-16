/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type Message } from '@molstar/ai-agent';
import { defaultDrive, drive } from '@molstar/ai-harness';
import { type ConformanceFixture, createSession, type Session, type SessionSetup } from '../session.js';
import { settle } from '../testing/clock.js';
import { system, textOf, user } from '../testing/messages.js';
import { expectWith, stepEvents, TURN, with_ } from './helpers.js';

/** Waits (in real turns of the event loop) until `ok` holds. */
async function until(ok: () => boolean | Promise<boolean>, rounds = 200): Promise<void> {
  for (let i = 0; i < rounds; i += 1) {
    if (await ok()) return;
    await settle();
  }
  throw new Error('until: the condition never held');
}

/** Steps (conformance.md): E1–E6. */
export function stepCases(fixture: ConformanceFixture): void {
  describe('Steps', () => {
    // E1. A step ends exactly one way, its commit matches that ending, and nothing is decided after the commit.

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

    it('E2 label-not-counted: a label never counts in a context size (it is never sent)', async () => {
      const sizes: number[] = [];
      const drive = with_({
        boundary: ({ payload }) => {
          sizes.push(payload.size.inputTokens);
          return 'pass';
        },
      });
      const prompt = user('go');
      const labelled: Message = { ...prompt, label: { kind: 'session:request', data: { note: 'x'.repeat(40_000) } } };
      for (const message of [prompt, labelled]) {
        const s = await createSession(fixture, { model: [{ text: 'ok' }], drive });
        await s.send([message]);
        expect((await s.main())[0]).toEqual(message);
      }
      expect(sizes).toHaveLength(2);
      expect(sizes[1]).toBe(sizes[0]);
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

    it('E4 retry-after: the wait a provider asks for is in the failure, its evidence, and the default retry', async () => {
      const asked: (number | undefined)[] = [];
      const s = await createSession(fixture, {
        model: [{ fail: { transient: true, retryAfterMs: 7_000 } }, { text: 'ok' }],
        drive: drive(
          {
            failure: ({ payload }) => {
              asked.push(payload.failure.retryAfterMs);
              return 'pass';
            },
          },
          defaultDrive({ retry: { max: 3, backoffMs: () => 10 } }),
        ),
      });
      const started = s.clock.now();
      await s.send('go');
      expect(asked).toEqual([7_000]);
      const [call] = stepEvents(await s.events()).flatMap((e) => e.calls.filter((c) => c.kind === 'model'));
      expect(call?.attempts[0]?.failure?.retryAfterMs).toBe(7_000);
      // The default drive waits what the provider asked, not its own shorter back-off.
      const second = call?.attempts[1];
      expect((second?.startedAt ?? 0) - started).toBeGreaterThanOrEqual(7_000);
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
  });
}
