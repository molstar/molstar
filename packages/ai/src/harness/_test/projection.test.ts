/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type Message } from '../../messages/index.js';
import { createSession } from '../../testing/session.js';
import { system, textOf, user } from '../../testing/messages.js';
import { stepEvents, with_ } from './helpers.js';

/** A small window, so a long prompt does not fit (harness.md §11, Admission). */
const SMALL_WINDOW = { contextWindowTokens: 10_000, maxOutputTokens: 1_000, defaultMaxOutputTokens: 1_000 };
const LONG = user('x'.repeat(60_000));

const toolIds = (tools: readonly { readonly tool: { readonly id: string } }[]) => tools.map((tool) => tool.tool.id);

/** The serialized messages up to, not including, the closing bracket: what a longer list must start with. */
const bytesOf = (messages: readonly Message[]) => JSON.stringify(messages).slice(0, -1);

/** C1–C2 (Projection). */
describe('Projection', () => {
  // C1 (harness.md §11)
  it('C1 same-input-same-digest: one head projected twice (in two harnesses) gives one digest', async () => {
    const [a, b] = [createSession({ model: [{ text: 'one' }] }), createSession({ model: [{ text: 'two' }] })];
    await a.send('go');
    await b.send('go');
    const [first, second] = [a.model.invocations[0], b.model.invocations[0]];
    expect(second?.digest).toBe(first?.digest);
    expect(second?.payload).toEqual(first?.payload);
    expect(b.model.projections[0]?.messages).toEqual(a.model.projections[0]?.messages);
  });

  it('C1 append-prefix: after a model-and-tool round in one turn, the earlier payload is a byte prefix of the new one', async () => {
    const s = createSession({
      tools: { look: () => ({ text: 'seen' }) },
      model: [{ calls: [{ tool: 'look', arguments: {} }] }, { text: 'done' }],
    });
    await s.send('go');
    const [first, second] = s.model.projections;
    expect(second?.messages.length).toBeGreaterThan(first?.messages.length ?? 0);
    expect(JSON.stringify(second?.messages).startsWith(bytesOf(first?.messages ?? []))).toBe(true);
    expect(toolIds(second?.tools ?? []).slice(0, toolIds(first?.tools ?? []).length)).toEqual(
      toolIds(first?.tools ?? []),
    );
  });

  // C2 (harness.md §11, Admission; §7.1)
  it('C2 over-window: failure with capacity opens before any attempt; compact applies the edit and the step projects again', async () => {
    const codes: string[] = [];
    const s = createSession({
      modelOptions: { capabilities: SMALL_WINDOW },
      model: [{ text: 'ok' }],
      drive: with_({
        failure: ({ payload }) => {
          codes.push(payload.failure.code);
          return { compact: { entries: [user('short')] } };
        },
      }),
    });
    await s.send([LONG]);
    expect(codes).toEqual(['capacity']);
    expect(s.model.projections).toHaveLength(2);
    expect(s.model.invocations).toHaveLength(1);
    expect(s.model.projections[1]?.messages).toEqual([user('short')]);
    expect(stepEvents(await s.events()).some((e) => e.compaction?.some((c) => c.pointer === 'main'))).toBe(true);
    const main = await s.main();
    expect(textOf(main[0] as Message)).toBe('short');
    expect((await s.turns())[0]?.status).toBe('completed');
  });

  it('C2 unrepresentable: a role the adapter cannot place yields unrepresentable, never a plain error (defect 15c)', async () => {
    const codes: string[] = [];
    const s = createSession({
      modelOptions: {
        unrepresentable: (input) => (input.messages.some((m) => m.role === 'system') ? 'no system role' : undefined),
      },
      model: [{ text: 'ok' }],
      drive: with_({
        failure: ({ payload }) => {
          codes.push(payload.failure.code);
          return 'pass';
        },
      }),
    });
    await s.send([system('be terse'), user('go')]);
    expect(codes).toEqual(['unrepresentable']);
    expect(s.model.invocations).toHaveLength(0);
    expect((await s.turns())[0]?.status).toBe('failed');
  });

  it('C2 unknown-tool-call: a call to a tool not in the request keeps its name and opens repair with unknown-tool (defect 15a)', async () => {
    const seen: { reason: string; toolId: string }[] = [];
    const s = createSession({
      model: [{ calls: [{ tool: 'ghost', arguments: {} }] }, { text: 'ok' }],
      drive: with_({
        repair: ({ payload }) => {
          seen.push({ reason: payload.reason, toolId: payload.call.toolId });
          return { report: { parts: [{ kind: 'text', text: 'no such tool' }] } };
        },
      }),
    });
    await s.send('go');
    expect(seen).toEqual([{ reason: 'unknown-tool', toolId: 'ghost' }]);
    const main = await s.main();
    const parts = main.flatMap((m): Message['parts'][number][] => [...m.parts]);
    const call = parts.find((p) => p.kind === 'tool-call');
    expect(call?.kind === 'tool-call' ? call.tool.id : undefined).toBe('ghost');
    const result = parts.find((p) => p.kind === 'tool-result');
    expect(result?.kind === 'tool-result' ? result.status : undefined).toBe('failed');
    expect((await s.turns())[0]?.status).toBe('completed');
  });
});
