/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message } from '@molstar/ai-agent';
import { type CallRecord, HarnessError, MAIN, type RunResult, type RunSpec } from '@molstar/ai-harness';
import { createSession, type Session } from '../session.js';
import { textOf, toolDefinition, user } from '../testing/messages.js';
import { stepEvents, until } from './helpers.js';

/** Handlers run inside the session they are made for; this is their way back to it. */
interface Later {
  s?: Session;
}

const nap = (later: Later, ms: number, signal: AbortSignal): Promise<void> =>
  later.s?.clock.sleep(ms, signal) ?? Promise.resolve();

/** A nested run over the given tools, all preauthorized, on the scripted model. */
const nestedRun = (...tools: string[]): RunSpec => ({
  input: [user('work on this')],
  tools: { enabled: tools, preauthorized: tools, definitions: tools.map((id) => toolDefinition(id)) },
  model: { target: { provider: 'scripted', model: 'scripted-1' }, parameters: {} },
  limits: {},
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
  const [pointer] = await s.harness.pointers([scope]);
  if (pointer === undefined) return [];
  const [head] = await s.harness.heads([pointer.head]);
  return head === undefined ? [] : await s.harness.entries(head.entries);
}

const nestedCalls = (calls: readonly CallRecord[]) => calls.filter((call) => call.parentCallId !== undefined);

/** N1–N3 (Nested runs). */
export function nestedCases(): void {
  describe('Nested runs', () => {
    it('N1 nested-records: a settled parent step holds the nested calls, with parentCallId (harness.md §9)', async () => {
      const s = createSession({
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
      const s = createSession({
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

    it('N2 pointer-busy: a nested run naming main fails pointer-busy (harness.md §8)', async () => {
      const outcomes: (string | undefined)[] = [];
      const s = createSession({
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

    it("N2 tool-journal: a handler's append moves its pointer when the step commits, not before (harness.md §8)", async () => {
      const later: Later = {};
      const reads: { before: number; after: number }[] = [];
      let reached = false;
      const s = createSession({
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
      const { handle } = await s.send('go', { wait: false });
      await until(() => reached);
      // The handler waits: its step is open, and the pointer has not moved.
      expect(await messagesOf(s, 'notes')).toEqual([]);
      await s.run(handle);
      // Within the step the handler read its own append.
      expect(reads).toEqual([{ before: 0, after: 1 }]);
      expect((await messagesOf(s, 'notes')).map(textOf)).toEqual(['remember']);
      const step = stepEvents(await s.events()).find((event) => event.stepKind === 'tools');
      expect(step?.moves.map((move) => move.pointer).sort()).toEqual([MAIN, 'notes']);
    });

    it("N3 nested-usage: the running budget excludes the nested 10k while the turn's usage totals include it", async () => {
      const later: Later = {};
      const seen: { consumed?: number | undefined; total?: number | undefined } = {};
      const s = createSession({
        configuration: { limits: { uncachedInput: 5_000 } },
        tools: {
          investigate: async (context) => ({ text: (await context.run(nestedRun())).status }),
          peek: async () => {
            const turn = (await later.s?.harness.execution())?.turn;
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
