/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message } from '@molstar/ai-agent';
import {
  type CallRecord,
  type Decision,
  HarnessError,
  type LaunchKey,
  MAIN,
  type RunnerExit,
  type RunResult,
  type RunSpec,
} from '@molstar/ai-harness';
import { type ConformanceFixture, createSession, type Session } from '../session.js';
import { settle } from '../testing/clock.js';
import { textOf, toolDefinition, user } from '../testing/messages.js';
import { stepEvents, with_ } from './helpers.js';

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

/** What a failed `run()` says, whether it rejects or returns a failed result (the contract says only that it fails). */
async function failureOf(attempt: Promise<RunResult>): Promise<string | undefined> {
  try {
    const result = await attempt;
    return result.status === 'failed' ? `failed: ${JSON.stringify(result.failure ?? null)}` : undefined;
  } catch (error) {
    return error instanceof HarnessError ? `${error.code}: ${error.message}` : String(error);
  }
}

/** The messages of a pointer's current head. */
async function messagesOf(s: Session, scope: string): Promise<Message[]> {
  const [pointer] = await s.client.pointers([scope]);
  if (pointer === undefined) return [];
  const [head] = await s.client.heads([pointer.head]);
  return head === undefined ? [] : await s.client.entries(head.entries);
}

const nestedCalls = (calls: readonly CallRecord[]) => calls.filter((call) => call.parentCallId !== undefined);

/** N1–N3 (Nested runs). */
export function nestedCases(fixture: ConformanceFixture): void {
  describe('Nested runs', () => {
    it('N1 nested-crash: a crash during a nested run writes nothing of it; the re-run runs it whole (harness.md §9)', async () => {
      const later: Later = {};
      const s = await createSession(fixture, {
        host: { leaseMs: 30_000 },
        declare: ['outer'],
        tools: {
          outer: async (context) => ({ text: (await context.run(nestedRun('inner'))).status }),
          inner: () => ({ text: 'inner done' }),
        },
        model: [
          { calls: [{ tool: 'outer', arguments: {} }] },
          { calls: [{ tool: 'inner', arguments: {} }] },
          { text: 'late', waitMs: 600_000 },
          { calls: [{ tool: 'inner', arguments: {} }] },
          { text: 'nested done' },
          { text: 'final' },
        ],
      });
      later.s = s;
      const { launch } = await s.send('go', { run: false });
      await crash(s, launch as string, () => s.model.invocations.length === 3);
      expect(s.tools.callsOf('inner')).toHaveLength(1);
      const before = stepEvents(await s.events()).flatMap((event) => event.calls);
      expect(nestedCalls(before)).toHaveLength(0);
      expect((await s.main()).some((message) => message.role === 'assistant')).toBe(false);
      await recover(s);
      expect(s.tools.callsOf('inner')).toHaveLength(2);
      const nested = nestedCalls(stepEvents(await s.events()).flatMap((event) => event.calls));
      expect(nested.filter((call) => call.kind === 'model')).toHaveLength(2);
      expect(nested.filter((call) => call.kind === 'tool')).toHaveLength(1);
      expect((await s.turns())[0]?.status).toBe('completed');
    });

    it('N1 nested-records: a settled parent step holds the nested calls, with parentCallId (harness.md §9)', async () => {
      const s = await createSession(fixture, {
        declare: ['outer'],
        tools: {
          outer: async (context) => ({ text: (await context.run(nestedRun('inner'))).status }),
          inner: () => ({ text: 'inner done' }),
        },
        model: [
          { calls: [{ tool: 'outer', arguments: {} }] },
          { calls: [{ tool: 'inner', arguments: {} }] },
          { text: 'nested done' },
          { text: 'final' },
        ],
      });
      await s.send('go');
      const outerId = s.tools.callsOf('outer')[0]?.callId;
      const events = stepEvents(await s.events());
      const step = events.find((event) => event.stepKind === 'tools');
      expect(step?.ending).toBe('settled');
      const outer = step?.calls.find((call) => call.id === outerId);
      expect(outer?.kind).toBe('tool');
      expect(outer?.parentCallId).toBeUndefined();
      const nested = step?.calls.filter((call) => call.parentCallId === outerId) ?? [];
      expect(nested.map((call) => call.kind).sort()).toEqual(['model', 'model', 'tool']);
      expect(step?.calls).toHaveLength(4);
      for (const call of nested) expect(call.attempts).toHaveLength(1);
      // No other step event carries nested calls: they belong to the parent's step only.
      expect(events.filter((event) => nestedCalls(event.calls).length > 0)).toEqual([step]);
    });

    it('N2 depth-exceeded: run() beyond depth 4 fails and the tool call reports it (harness.md §9)', async () => {
      const outcomes: (string | undefined)[] = [];
      const s = await createSession(fixture, {
        declare: ['dive'],
        tools: {
          dive: async (context) => {
            const failure = await failureOf(context.run(nestedRun('dive')));
            outcomes.push(failure);
            return failure === undefined ? { text: 'ok' } : { status: 'failed', text: `run refused: ${failure}` };
          },
        },
        model: [
          // The top level and the nested runs at depths 1 to 4 each call dive; the fifth run() is refused.
          ...Array.from({ length: 5 }, () => ({ calls: [{ tool: 'dive', arguments: {} }] })),
          { text: 'depth 4 done' },
          { text: 'depth 3 done' },
          { text: 'depth 2 done' },
          { text: 'depth 1 done' },
          { text: 'final' },
        ],
      });
      await s.send('go');
      expect(s.tools.callsOf('dive')).toHaveLength(5);
      expect(outcomes).toHaveLength(5);
      expect(outcomes[0]).toBeDefined();
      expect(outcomes.slice(1)).toEqual([undefined, undefined, undefined, undefined]);
      expect(s.model.remaining()).toBe(0);
      expect((await s.turns())[0]?.status).toBe('completed');
    });

    it('N2 pointer-busy: of two sibling nested runs naming one pointer the second fails pointer-busy (harness.md §8, §9)', async () => {
      const outcomes: (string | undefined)[] = [];
      const s = await createSession(fixture, {
        declare: ['a', 'b'],
        tools: {
          a: async (context) => {
            const failure = await failureOf(context.run({ ...nestedRun(), pointer: 'shared' }));
            outcomes.push(failure);
            return { text: failure ?? 'ran' };
          },
          b: async (context) => {
            const failure = await failureOf(context.run({ ...nestedRun(), pointer: 'shared' }));
            outcomes.push(failure);
            return { text: failure ?? 'ran' };
          },
        },
        model: [
          {
            calls: [
              { tool: 'a', arguments: {} },
              { tool: 'b', arguments: {} },
            ],
          },
          { text: 'nested done' },
          { text: 'final' },
        ],
        drive: with_({ schedule: parallel }),
      });
      await s.send('go');
      expect(outcomes).toHaveLength(2);
      expect(outcomes.filter((outcome) => outcome === undefined)).toHaveLength(1);
      expect(outcomes.filter((outcome) => outcome?.includes('pointer-busy'))).toHaveLength(1);
      expect(await messagesOf(s, 'shared')).toHaveLength(2);
      expect(s.model.remaining()).toBe(0);
    });

    it('N2 pointer-busy: a nested run naming main fails pointer-busy (harness.md §8)', async () => {
      const outcomes: (string | undefined)[] = [];
      const s = await createSession(fixture, {
        declare: ['m'],
        tools: {
          m: async (context) => {
            const failure = await failureOf(context.run({ ...nestedRun(), pointer: MAIN }));
            outcomes.push(failure);
            return { text: failure ?? 'ran' };
          },
        },
        model: [{ calls: [{ tool: 'm', arguments: {} }] }, { text: 'final' }],
      });
      await s.send('go');
      expect(outcomes).toHaveLength(1);
      expect(outcomes[0]).toContain('pointer-busy');
      expect(s.model.remaining()).toBe(0);
    });

    it('N2 fork: a nested run started from a head moves no pointer (harness.md §9)', async () => {
      const s = await createSession(fixture, {
        declare: ['fork'],
        tools: {
          fork: async (context) => {
            const { head } = await context.journal.read(MAIN);
            if (head === undefined) throw new Error('main has no head');
            const result = await context.run({ ...nestedRun(), input: [user('fork input')], from: head });
            return { text: result.status };
          },
        },
        model: [{ calls: [{ tool: 'fork', arguments: {} }] }, { text: 'forked' }, { text: 'final' }],
      });
      await s.send('go');
      const step = stepEvents(await s.events()).find((event) => event.stepKind === 'tools');
      expect(step?.moves.map((move) => move.pointer)).toEqual([MAIN]);
      expect((await s.main()).map(textOf)).not.toContain('fork input');
      // The fork continued from main's head: its model saw main's prompt before its own input.
      const texts = s.model.projections[1]?.messages.map(textOf) ?? [];
      expect(texts).toContain('go');
      expect(texts).toContain('fork input');
      expect(texts.indexOf('go')).toBeLessThan(texts.indexOf('fork input'));
    });

    it("N2 tool-journal: a handler's append moves its pointer when the step commits, not before (harness.md §8)", async () => {
      const later: Later = {};
      const reads: { before: number; after: number }[] = [];
      let reached = false;
      const s = await createSession(fixture, {
        host: { leaseMs: 30_000 },
        tools: {
          note: async (context) => {
            const before = (await context.journal.read('notes')).messages.length;
            context.journal.append('notes', [user('remember')]);
            const after = (await context.journal.read('notes')).messages.map(textOf);
            reads.push({ before, after: after.filter((text) => text === 'remember').length });
            if (reads.length === 1) {
              reached = true;
              await nap(later, 600_000, context.signal);
            }
            return { text: 'noted' };
          },
        },
        model: [{ calls: [{ tool: 'note', arguments: {} }] }, { text: 'final' }],
      });
      later.s = s;
      const { launch } = await s.send('go', { run: false });
      await crash(s, launch as string, () => reached);
      expect(await messagesOf(s, 'notes')).toEqual([]);
      await recover(s);
      // Within the step the handler read its own append; the crashed attempt moved nothing.
      expect(reads).toEqual([
        { before: 0, after: 1 },
        { before: 0, after: 1 },
      ]);
      expect((await messagesOf(s, 'notes')).map(textOf)).toEqual(['remember']);
      const step = stepEvents(await s.events()).find((event) => event.stepKind === 'tools');
      expect(step?.moves.map((move) => move.pointer).sort()).toEqual([MAIN, 'notes']);
    });

    it("N2 main-read-only: a handler's append to main fails pointer-busy; read works (harness.md §8)", async () => {
      let code: string | undefined;
      let readable = 0;
      const s = await createSession(fixture, {
        tools: {
          sneak: async (context) => {
            try {
              context.journal.append(MAIN, [user('sneaky')]);
            } catch (error) {
              code = error instanceof HarnessError ? error.code : String(error);
            }
            readable = (await context.journal.read(MAIN)).messages.length;
            return { text: 'tried' };
          },
        },
        model: [{ calls: [{ tool: 'sneak', arguments: {} }] }, { text: 'final' }],
      });
      await s.send('go');
      expect(code).toBe('pointer-busy');
      expect(readable).toBeGreaterThan(0);
      expect((await s.main()).map(textOf)).not.toContain('sneaky');
    });

    it("N2 hold-per-attempt: a failed attempt's append is dropped; its retry appends again and only that commits (harness.md §8)", async () => {
      const errors: string[] = [];
      const s = await createSession(fixture, {
        tools: {
          log: (context) => {
            try {
              context.journal.append('log', [user(`attempt ${context.attempt}`)]);
            } catch (error) {
              errors.push(error instanceof HarnessError ? error.code : String(error));
            }
            return context.attempt === 1 ? { status: 'failed', text: 'flaky', transient: true } : { text: 'logged' };
          },
        },
        model: [{ calls: [{ tool: 'log', arguments: {} }] }, { text: 'final' }],
      });
      await s.send('go');
      expect(s.tools.callsOf('log').map((call) => call.attempt)).toEqual([1, 2]);
      expect(errors).toEqual([]);
      expect((await messagesOf(s, 'log')).map(textOf)).toEqual(['attempt 2']);
    });

    it("N2 pointer-continues: a nested run's pointer, used on a later turn, continues from its head (harness.md §9)", async () => {
      let uses = 0;
      const s = await createSession(fixture, {
        declare: ['helper'],
        tools: {
          helper: async (context) => {
            uses += 1;
            const input = [user(uses === 1 ? 'first' : 'second')];
            return { text: (await context.run({ ...nestedRun(), input, pointer: 'helper' })).status };
          },
        },
        model: [
          { calls: [{ tool: 'helper', arguments: {} }] },
          { text: 'one' },
          { text: 'turn 1 done' },
          { calls: [{ tool: 'helper', arguments: {} }] },
          { text: 'two' },
          { text: 'turn 2 done' },
        ],
      });
      await s.send('go');
      expect((await messagesOf(s, 'helper')).map(textOf)).toEqual(['first', 'one']);
      await s.send('again');
      expect(await s.turns()).toHaveLength(2);
      expect((await messagesOf(s, 'helper')).map(textOf)).toEqual(['first', 'one', 'second', 'two']);
      // The second turn's nested model saw the first turn's exchange before the new input.
      const seen = s.model.projections.find((projection) => projection.messages.map(textOf).includes('second'));
      expect(seen?.messages.map(textOf).filter((text) => ['first', 'one', 'second'].includes(text))).toEqual([
        'first',
        'one',
        'second',
      ]);
    });

    it("N3 nested-usage: a nested run's usage is reported, never charged to the parent", async () => {
      const s = await createSession(fixture, {
        configuration: { limits: { uncachedInput: 5_000 } },
        tools: {
          investigate: async (context) => {
            const result = await context.run({
              input: [user('look closely')],
              tools: { enabled: [], preauthorized: [], definitions: [] },
              model: { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} },
              limits: {},
            });
            return { text: result.status };
          },
        },
        model: [
          { calls: [{ tool: 'investigate', arguments: {} }] },
          { text: 'found it', usage: { inputTokens: 10_000, outputTokens: 10 } },
          { text: 'final' },
        ],
      });
      await s.send('go');
      const [turn] = await s.turns();
      expect(turn?.status).toBe('completed');
      expect(turn?.usage.uncachedInput).toBeGreaterThanOrEqual(10_000);
    });

    it("N3 nested-usage: the running budget excludes the nested 10k while the turn's usage totals include it", async () => {
      const later: Later = {};
      const seen: { consumed?: number | undefined; total?: number | undefined } = {};
      const s = await createSession(fixture, {
        configuration: { limits: { uncachedInput: 5_000 } },
        tools: {
          investigate: async (context) => ({ text: (await context.run(nestedRun())).status }),
          peek: async () => {
            const turn = (await later.s?.client.execution())?.turn;
            seen.consumed = turn?.budget.consumed.uncachedInput;
            seen.total = turn?.usage.uncachedInput;
            return { text: 'peeked' };
          },
        },
        model: [
          { calls: [{ tool: 'investigate', arguments: {} }] },
          { text: 'found it', usage: { inputTokens: 10_000, outputTokens: 10 } },
          { calls: [{ tool: 'peek', arguments: {} }] },
          { text: 'final' },
        ],
      });
      later.s = s;
      await s.send('go');
      // The first round (nested run included) is committed by the time the second round's tool looks.
      expect(seen.consumed).toBeDefined();
      expect(seen.consumed).toBeLessThan(10_000);
      expect(seen.total).toBeGreaterThanOrEqual(10_000);
      const step = stepEvents(await s.events()).find((event) => event.stepKind === 'tools');
      const nestedInput = nestedCalls(step?.calls ?? []).reduce(
        (sum, call) => sum + call.attempts.reduce((n, attempt) => n + attempt.usage.uncachedInput, 0),
        0,
      );
      expect(nestedInput).toBeGreaterThanOrEqual(10_000);
    });
  });
}
