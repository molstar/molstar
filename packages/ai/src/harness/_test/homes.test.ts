/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message, ToolDefinition } from '../../messages/index.js';
import { createSession } from '../../testing/session.js';
import { assistant, declareTools, textOf, toolDefinition } from '../../testing/messages.js';
import { with_ } from './helpers.js';

/** H1, all negative: each concern stays in its home (harness.md §1). */
describe('Homes', () => {
  it('H1 contribution-like-answer: an assistant-text contribution does not end the turn (defect 2)', async () => {
    const s = createSession({
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
    const s = createSession({
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

  it('H1 redeclared: a later redeclaration of a tool does not change the definition a planned call runs under (defect 7)', async () => {
    const redeclared = toolDefinition('probe', { depth: { type: 'string' } });
    const seen: ToolDefinition[] = [];
    const s = createSession({
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
});
