/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { digest, hashMessage, type Json, type Message, type ModelAdapter } from '../../messages/index.js';
import { createHarness, defaultDrive, emptyHead, MAIN, type TurnConfiguration } from '../index.js';

// The harness's state in memory (harness.md §2–§4): what the trace and the Inspector read of it.

const user = (text: string): Message => ({ role: 'user', parts: [{ kind: 'text', text }] });

/** A model that answers each request with the next of `replies`. */
function replying(replies: string[]): ModelAdapter {
  return {
    capabilitiesFor: async () => ({
      contextWindowTokens: 100_000,
      maxOutputTokens: 1_000,
      outputLimitIncludesReasoning: false,
      appendsTools: true,
      cacheScope: true,
    }),
    async project(input) {
      const payload = input.messages.map(({ role, parts }) => ({ role, parts })) as unknown as Json;
      return { payload, digest: digest(payload), adapterRevision: 'test-1', inputTokens: 10 };
    },
    async *invoke() {
      const text = replies.shift() ?? 'nothing left';
      yield {
        kind: 'result',
        message: { role: 'assistant', parts: [{ kind: 'text', text }] },
        finish: 'complete',
        usage: { inputTokens: 10, outputTokens: 2 },
      };
    },
  };
}

const configuration: TurnConfiguration = {
  model: { target: { provider: 'test', model: 'test-1' }, parameters: {} },
  tools: { enabled: [], preauthorized: [] },
  limits: {},
};

function harness(replies: string[]) {
  const adapter = replying(replies);
  return createHarness({
    id: 'h-state',
    models: { resolve: async () => adapter },
    tools: { resolve: () => undefined },
    drive: defaultDrive(),
    configure: () => configuration,
  });
}

describe('the harness state', () => {
  it('starts empty: main is the empty head, no events, no turns', async () => {
    const h = harness([]);
    expect(await h.execution()).toEqual({ harnessId: 'h-state', seq: 0, main: emptyHead().id });
    expect(await h.heads([emptyHead().id])).toEqual([emptyHead()]);
    expect(await h.events({})).toEqual([]);
    expect(await h.turns({})).toEqual([]);
    expect(await h.pointers([MAIN, 'notes'])).toEqual([{ scope: MAIN, head: emptyHead().id }, undefined]);
  });

  it('keeps an identical message once, and heads as lists of entry ids', async () => {
    const h = harness(['one', 'two']);
    await h.turn({ messages: [user('same')] }).done;
    await h.turn({ messages: [user('same')] }).done;
    const { main } = await h.execution();
    const [head] = await h.heads([main]);
    const id = hashMessage(user('same'));
    expect(head?.entries.filter((entry) => entry === id)).toHaveLength(2);
    expect(new Set(head?.entries).size).toBe(3);
    expect((await h.entries([id]))[0]).toEqual(user('same'));
    await expect(h.entries(['no-such-entry'])).rejects.toMatchObject({ code: 'invalid' });
  });

  it('pages the event log by seq, reads steps by id and turns newest first', async () => {
    const h = harness(['one', 'two']);
    const first = h.turn({ messages: [user('a')] });
    await first.done;
    const second = h.turn({ messages: [user('b')] });
    await second.done;
    const events = await h.events({ limit: 100 });
    expect(events.map((event) => event.seq)).toEqual(events.map((_, i) => i + 1));
    expect(await h.events({ after: 2, limit: 1 })).toEqual([events[2]]);
    const turns = await h.turns({});
    expect(turns.map((turn) => [turn.id, turn.status])).toEqual([
      [second.turnId, 'completed'],
      [first.turnId, 'completed'],
    ]);
    // Before the second turn's last event: the first turn only.
    expect(await h.turns({ before: events.at(-1)?.seq, limit: 1 })).toEqual([turns[1]]);
    const [ending] = await h.steps([turns[0]?.lastStep ?? '']);
    expect(ending?.end).toEqual({ status: 'completed' });
    const [request] = await h.blobs([ending?.calls[0]?.request ?? '']);
    expect(request).toMatchObject({ adapterRevision: 'test-1' });
  });

  it('watch gives the current seq first, then each change', async () => {
    const h = harness(['one']);
    const seen: number[] = [];
    const stop = h.watch().subscribe((seq) => seen.push(seq));
    await h.turn({ messages: [user('a')] }).done;
    stop();
    expect(seen).toEqual([0, 1, 2]);
  });
});
