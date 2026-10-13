/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message, ToolDefinition } from '@molstar/ai-agent';
import type { LastCompaction } from '@molstar/ai-harness';
import { describe, expect, it } from '@jest/globals';
import { type ConformanceFixture, createSession } from '../session.js';
import { settle } from '../testing/clock.js';
import { assistant, declareTools, textOf, toolDefinition, user } from '../testing/messages.js';
import { stepEvents, TURN, with_ } from './helpers.js';

/** Waits (in real turns of the event loop) until `ok` holds. */
async function until(ok: () => boolean | Promise<boolean>, rounds = 200): Promise<void> {
  for (let i = 0; i < rounds; i += 1) {
    if (await ok()) return;
    await settle();
  }
  throw new Error('until: the condition never held');
}

/** H1, all negative: each concern stays in its home (harness.md §1). */
export function homeCases(fixture: ConformanceFixture): void {
  describe('Homes', () => {
    it('H1 contribution-like-answer: an assistant-text contribution does not end the turn (defect 2)', async () => {
      const s = await createSession(fixture, {
        tools: { note: () => ({ text: 'noted', contributions: [assistant('All done.')] }) },
        model: [{ calls: [{ tool: 'note', arguments: {} }] }, { text: 'final' }],
      });
      await s.send('go');
      expect(s.model.invocations).toHaveLength(2);
      const main = await s.main();
      const last = main.at(-1);
      expect(last?.role).toBe('assistant');
      expect(last === undefined ? '' : textOf(last)).toBe('final');
    });

    it('H1 contribution-with-calls: a contribution holding an unanswered call fails invalid-contribution and plans nothing (defect 2)', async () => {
      const smuggled: Message = {
        role: 'assistant',
        parts: [{ kind: 'tool-call', requestId: 'smuggled', tool: { id: 'note' }, arguments: {} }],
      };
      const failures: string[] = [];
      const s = await createSession(fixture, {
        tools: { note: () => ({ text: 'noted', contributions: [smuggled] }) },
        model: [{ calls: [{ tool: 'note', arguments: {} }] }, { text: 'final' }],
        drive: with_({
          failure: ({ payload }) => {
            failures.push(payload.failure.code);
            return { report: {} };
          },
        }),
      });
      await s.send('go');
      // harness.md §8: the attempt fails invalid-contribution (harness.md §16: an attempt failure opens `failure`).
      expect(failures).toEqual(['invalid-contribution']);
      // Nothing is planned from the contribution: the handler ran once, and the model was asked once more.
      expect(s.tools.callsOf('note')).toHaveLength(1);
      expect(s.model.invocations).toHaveLength(2);
      const main = await s.main();
      expect(JSON.stringify(main)).not.toContain('smuggled');
      const results = main.flatMap((m) => m.parts.flatMap((p) => (p.kind === 'tool-result' ? [p] : [])));
      expect(results.map((r) => r.status)).toEqual(['failed']);
    });

    it('H1 edited-head: a compaction edit that removes the last assistant message does not change the named next step', async () => {
      let compacted = false;
      const s = await createSession(fixture, {
        model: [{ text: 'final' }, { text: 'must not be asked for' }],
        drive: with_({
          boundary: ({ payload }, ctx) => {
            if (!payload.completing || compacted) return 'proceed';
            compacted = true;
            // At the completing boundary the transient head ends with the assistant message just made.
            return { compact: { entries: ctx.head.entries.slice(0, -1) } };
          },
        }),
      });
      await s.send('go');
      expect(compacted).toBe(true);
      // The head now ends with a user message, which a head-derived rule would answer with a model
      // step. The last step event named `end: completed` (harness.md §6: never read from a head).
      expect(s.model.invocations).toHaveLength(1);
      const [turn] = await s.turns();
      expect(turn?.status).toBe('completed');
      const last = stepEvents(await s.events()).at(-1);
      expect(last?.end?.status).toBe('completed');
      expect(last?.next).toBeUndefined();
    });

    it('H1 redeclared: a later redeclaration of a tool does not change the definition a planned call runs under (defect 7)', async () => {
      const redeclared = toolDefinition('probe', { depth: { type: 'string' } });
      const seen: ToolDefinition[] = [];
      const s = await createSession(fixture, {
        tools: {
          probe: (context) => {
            seen.push(context.definition);
            return { text: 'ok', contributions: [declareTools([redeclared])] };
          },
        },
        model: [
          { calls: [{ tool: 'probe', arguments: {} }] },
          { calls: [{ tool: 'probe', arguments: {} }] },
          { text: 'final' },
        ],
      });
      await s.send('go');
      // harness.md §8: a call binds to the last loaded definition that precedes its assistant message.
      expect(seen).toEqual([s.tools.definition('probe'), redeclared]);
    });

    it('H1 no-evidence-entries: a step with calls and a failed retry leaves no request, attempt or failure in main', async () => {
      const s = await createSession(fixture, {
        tools: { note: () => ({ text: 'noted' }) },
        model: [
          { fail: { message: 'evidence-boom', transient: true } },
          { calls: [{ tool: 'note', arguments: {} }] },
          { text: 'final' },
        ],
      });
      await s.send('go');
      const steps = stepEvents(await s.events());
      const modelCalls = steps.flatMap((e) => e.calls.filter((c) => c.kind === 'model'));
      // The evidence exists in the event log (harness.md §3) ...
      expect(modelCalls.some((c) => c.attempts.some((a) => a.outcome === 'failed'))).toBe(true);
      // ... and none of it is in the journal.
      const main = JSON.stringify(await s.main());
      const evidence = steps.flatMap((e) => e.calls.flatMap((c) => [c.id, c.request, ...c.attempts.map((a) => a.id)]));
      expect(evidence.length).toBeGreaterThan(0);
      for (const id of evidence) expect(main).not.toContain(id);
      expect(main).not.toContain('evidence-boom');
    });

    it('H1 last-compaction-after-stop: compact, Stop, a new turn; next-turn and boundary still report the compaction (defect 11)', async () => {
      const seen: { kind: string; last: LastCompaction | undefined }[] = [];
      let compacted = false;
      const s = await createSession(fixture, {
        tools: { note: () => ({ text: 'noted' }) },
        model: [
          { calls: [{ tool: 'note', arguments: {} }] },
          { text: 'slow', waitMs: 60_000 },
          { text: 'second turn' },
        ],
        drive: with_({
          boundary: ({ payload }, ctx) => {
            seen.push({ kind: 'boundary', last: payload.lastCompaction });
            if (compacted) return 'proceed';
            compacted = true;
            return { compact: { entries: [...ctx.head.entries.slice(0, 1), user('summary of the work so far')] } };
          },
          'next-turn': ({ payload }) => {
            seen.push({ kind: 'next-turn', last: payload.lastCompaction });
            return 'pass';
          },
        }),
      });
      const { launch } = await s.send('go', { run: false });
      const running = s.runner().run(launch as string);
      await until(() => s.model.invocations.length === 2);
      const turnId = (await s.client.execution()).turn?.id as string;
      await s.client.stop({ key: 'stop-1', turnId });
      await s.clock.runUntil(running);
      expect((await s.turns()).map((t) => t.status)).toContain('stopped');

      await s.send('again');
      const compaction = stepEvents(await s.events()).find((e) => (e.compaction?.length ?? 0) > 0);
      const size = compaction?.compaction?.[0]?.sizeAfter;
      expect(size).toBeDefined();
      const expected = { seq: compaction?.seq, sizeAfter: size };
      // The new turn's first next-turn and its completing boundary both report it (harness.md §3, §7).
      expect(seen.filter((x) => x.kind === 'next-turn').at(-1)?.last).toEqual(expected);
      expect(seen.filter((x) => x.kind === 'boundary').at(-1)?.last).toEqual(expected);
    });

    it("H1 compaction-per-pointer: lastCompaction(main) is main's when a nested run also compacted its own pointer", async () => {
      let nestedCompacted = false;
      let mainCompacted = false;
      const seen: { completing: boolean; last: LastCompaction | undefined }[] = [];
      const s = await createSession(fixture, {
        tools: {
          helper: async (context) => {
            const result = await context.run({
              input: [user('help')],
              tools: { enabled: [], preauthorized: [], definitions: [] },
              model: TURN.model,
              pointer: 'helper',
              limits: {},
              drive: with_({
                boundary: () => {
                  if (nestedCompacted) return 'proceed';
                  nestedCompacted = true;
                  return { compact: { entries: [user('h'.repeat(4_000))] } };
                },
              }),
            });
            return { text: result.status };
          },
        },
        model: [{ calls: [{ tool: 'helper', arguments: {} }] }, { text: 'nested output' }, { text: 'final' }],
        drive: with_({
          boundary: ({ payload }, ctx) => {
            seen.push({ completing: payload.completing, last: payload.lastCompaction });
            if (payload.completing || mainCompacted) return 'proceed';
            mainCompacted = true;
            return { compact: { entries: [...ctx.head.entries.slice(0, 1), user('main summary')] } };
          },
        }),
      });
      await s.send('go');
      const step = stepEvents(await s.events()).find((e) => (e.compaction?.length ?? 0) > 1);
      expect(step?.compaction?.map((c) => c.pointer).sort()).toEqual(['helper', 'main']);
      const main = step?.compaction?.find((c) => c.pointer === 'main');
      const helper = step?.compaction?.find((c) => c.pointer === 'helper');
      expect(helper?.sizeAfter).not.toEqual(main?.sizeAfter);
      // After that step's commit, the completing boundary reports main's compaction, not the helper's.
      expect(seen.at(-1)).toEqual({ completing: true, last: { seq: step?.seq, sizeAfter: main?.sizeAfter } });
    });

    it.todo("H1 import-rule: Biome's (biome.json, noRestrictedImports); it checks source imports, not behavior");
  });
}
