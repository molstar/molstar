/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createMemoryStorage } from '@molstar/ai-components';
import { createToolkit } from '@molstar/ai-tools';
import { GUIDANCE_UPDATE, type StatusEvent, type ToolEvent, type TurnEvent } from '../index.js';
import { callOf, fixture, guidance, headOf, messagesOf, realClock, resultsIn, textOf, tool } from './helpers.js';

describe('a turn (tools-and-session.md §4)', () => {
  it('runs a tool call and answers in words; the hooks see the text, the calls, the status and the turn', async () => {
    const { session, api, model } = fixture([{ calls: [callOf('set_value', { by: '2' })] }, { text: 'Bumped it.' }]);
    const texts: string[] = [];
    const calls: ToolEvent[] = [];
    const statuses: StatusEvent['kind'][] = [];
    const turns: TurnEvent[] = [];
    const order: string[] = [];
    const r = await session.turn(
      { prompt: 'bump the value by two' },
      {
        onText: (t) => texts.push(t),
        onTool: (e) => {
          calls.push(e);
          order.push(`tool:${e.status}`);
        },
        onStatus: (s) => statuses.push(s.kind),
        onTurn: (e) => {
          turns.push(e);
          order.push('turn');
        },
      },
    );
    expect(r.status).toBe('completed');
    expect(r.text).toBe('Bumped it.');
    expect(texts.at(-1)).toBe('Bumped it.');
    expect(api.state.value).toBe(2);
    expect(r.tools.map((x) => [x.name, x.status, x.read])).toEqual([['set_value', 'succeeded', false]]);
    expect(r.tools[0]?.text).toBe('Bumped.\nChanged: value.');
    expect(calls.map((x) => x.status)).toEqual(['started', 'succeeded']);
    // Times come from the session's clock.
    const done = r.tools[0];
    expect(done?.startedAt).toBeGreaterThan(0);
    expect(done?.endedAt).toBeGreaterThanOrEqual(done?.startedAt ?? Number.POSITIVE_INFINITY);
    // The turn is announced once, before any call starts.
    expect(turns).toEqual([{ kind: 'begun', turnId: expect.any(String) }]);
    expect(order[0]).toBe('turn');
    expect(statuses).toEqual(expect.arrayContaining(['model', 'tools']));
    expect(r.usage?.modelCalls).toBe(2);
    expect(r.context).toBeGreaterThan(0);
    expect(model.remaining()).toBe(0);
    expect(session.busy).toBe(false);
  });

  it('marks a read, and sends a failure to the model, which corrects itself', async () => {
    const { session, api } = fixture([
      { calls: [callOf('get_state'), callOf('set_value', { by: 'many' })] },
      { calls: [callOf('set_value', { by: 3 })] },
      { text: 'Done.' },
    ]);
    const r = await session.turn({ prompt: 'bump it' });
    expect(r.tools.map((x) => [x.name, x.status, x.read])).toEqual([
      ['get_state', 'succeeded', true],
      ['set_value', 'failed', false],
      ['set_value', 'succeeded', false],
    ]);
    expect(r.tools[1]?.text).toBe('by: expected number (got "many").');
    expect(api.state.value).toBe(3);
  });

  it('sends one turn at a time: a second is refused while one runs', async () => {
    const { session } = fixture([{ text: 'Done.' }]);
    const first = session.turn({ prompt: 'one' });
    expect(session.busy).toBe(true);
    await expect(session.turn({ prompt: 'two' })).rejects.toThrow('A turn is already running.');
    expect((await first).status).toBe('completed');
    expect(session.busy).toBe(false);
  });

  it('a tool loaded by name is then called; a skill loads its tools with its guidance', async () => {
    const { session, api } = fixture([
      { calls: [callOf('disclose_tools', { tools: ['archive'] })] },
      { calls: [callOf('archive'), callOf('disclose_skills', { skills: ['maps'] })] },
      { calls: [callOf('load_map')] },
      { text: 'Done.' },
    ]);
    const r = await session.turn({ prompt: 'archive it and load a map' });
    expect(r.status).toBe('completed');
    expect(r.tools.map((x) => `${x.name}:${x.status}`)).toEqual([
      'disclose_tools:succeeded',
      'archive:succeeded',
      'disclose_skills:succeeded',
      'load_map:succeeded',
    ]);
    expect(r.tools[0]?.text).toBe('Loaded: archive.');
    expect(api.log.lines).toEqual(['load_map']);
    // The guidance of the skill is in the conversation, as its own entry.
    expect((await headOf(session)).some((m) => textOf(m).startsWith('Skill maps: load a map.'))).toBe(true);
  });

  it('stops a running turn', async () => {
    const { session } = fixture([{ text: 'Never said.', waitMs: 10_000 }], {}, { clock: realClock });
    const running = session.turn({ prompt: 'wait' });
    await new Promise((resolve) => setTimeout(resolve, 100));
    await session.stop();
    const r = await running;
    expect(r.status).toBe('stopped');
    expect(session.busy).toBe(false);
  });

  it('repeats a read briefly while the state is unchanged, in full once it changes; a long streak is nudged', async () => {
    const { session, api } = fixture(
      [
        { calls: [callOf('get_state')] },
        { calls: [callOf('get_state')] },
        { calls: [callOf('set_value')] },
        { calls: [callOf('get_state')] },
        ...[1, 2, 3, 4, 5, 6].map((look) => ({ calls: [callOf('get_state', { look })] })),
        { text: 'Done.' },
      ],
      { version: () => api.state.version },
      {},
    );
    const r = await session.turn({ prompt: 'look' });
    const reads = r.tools.filter((t) => t.name === 'get_state').map((t) => t.text ?? '');
    expect(reads.map((t) => t.startsWith('Unchanged'))).toEqual([
      false,
      true,
      false,
      false,
      false,
      false,
      false,
      false,
      false,
    ]);
    // From the sixth read in a row since the change on, each result says to act.
    expect(reads.map((t) => t.includes('reads in a row'))).toEqual([
      false,
      false,
      false,
      false,
      false,
      false,
      false,
      true,
      true,
    ]);
  });
});

describe('what the host hooks into', () => {
  it('saves after a write before its result is kept, not after a read; a failed save is not retried', async () => {
    const saved: string[] = [];
    const f = fixture([{ calls: [callOf('get_state'), callOf('set_value')] }, { text: 'Done.' }], {
      afterWrite: async ({ toolId }) => void saved.push(toolId),
    });
    await f.session.turn({ prompt: 'read and bump' });
    expect(saved).toEqual(['set_value']);

    const broken = fixture([{ calls: [callOf('set_value')] }, { text: 'Not saved.' }], {
      afterWrite: async () => {
        throw new Error('disk full');
      },
    });
    const r = await broken.session.turn({ prompt: 'bump' });
    // The change was made, so running it again could repeat it: the failure goes to the model.
    expect(broken.api.state.value).toBe(1);
    expect(r.tools.map((t) => [t.name, t.status])).toEqual([['set_value', 'failed']]);
    expect(r.tools[0]?.text).toBe('The change could not be saved: disk full');
  });

  it('prices the model calls with the host’s function', async () => {
    const f = fixture([{ text: 'Done.' }], { price: (_configuration, usage) => usage.modelCalls * 0.5 });
    await f.session.turn({ prompt: 'hi' });
    const client = await f.session.client();
    expect((await client.turns({ limit: 1 }))[0]?.usage.cost).toBe(0.5);
  });

  it('knows the default model’s window: the host’s, else what its adapter reported', async () => {
    const given = fixture([{ text: 'Done.' }], { contextWindow: 500_000 });
    expect(given.session.contextWindow()).toBe(500_000);
    const reported = fixture([{ text: 'Done.' }], {}, { capabilities: { contextWindowTokens: 300_000 } });
    expect(reported.session.contextWindow()).toBe(200_000);
    await reported.session.turn({ prompt: 'hi' });
    expect(reported.session.contextWindow()).toBe(300_000);
  });

  it('heads the notes with the host’s words', async () => {
    const f = fixture([{ text: 'Done.' }], { notesHeading: 'The person changed:' });
    await f.session.turn({ prompt: 'hi', notes: ['the color'] });
    expect(textOf(messagesOf(f.model, 0).at(-1) as never)).toBe('The person changed:\n- the color\n\nRequest: hi');
  });

  it('shows the person a tool’s own calls, not those of a helper it runs', async () => {
    const inner = tool({
      id: 'inner',
      description: 'The helper’s tool.',
      input: { kind: 'object', fields: {} },
      modules: ['state'],
      run: ({ api }) => {
        api.state.bump(10);
        return { ok: true, did: 'Inner done.' };
      },
    });
    const helper = tool({
      id: 'helper',
      description: 'Runs a helper.',
      input: { kind: 'object', fields: {} },
      modules: [],
      run: async ({ nest }) => {
        const r = await nest?.({
          messages: [{ role: 'user', parts: [{ kind: 'text', text: 'Help.' }] }],
          tools: [inner],
          api: f.api,
          model: { provider: 'scripted', model: 'scripted-1' },
          limits: {},
        });
        return { ok: true, did: `Helper ${r?.status}.` };
      },
    });
    const toolkit = createToolkit<typeof f.api>({ tools: [helper], core: ['helper'], groups: [] });
    const f = fixture(
      [{ calls: [callOf('helper')] }, { calls: [callOf('inner')] }, { text: 'Helped.' }, { text: 'All done.' }],
      { toolkit },
    );
    const seen: string[] = [];
    const r = await f.session.turn({ prompt: 'use the helper' }, { onTool: (e) => seen.push(`${e.name}:${e.status}`) });
    expect(r.status).toBe('completed');
    expect(f.api.state.value).toBe(10);
    expect(seen).toEqual(['helper:started', 'helper:succeeded']);
    expect(r.tools.map((t) => t.text)).toEqual(['Helper completed.']);
  });
});

describe('the messages (tools-and-session.md §4)', () => {
  it('the first turn sends the guidance, the tools (the core loaded, the rest declared), then the request', async () => {
    const { session, model, toolkit } = fixture([{ text: 'Hi.' }]);
    await session.turn({
      prompt: 'hello',
      state: 'Counter at 0.',
      focus: 'the display',
      notes: ['The person reset the counter.'],
      preamble: 'Memory: likes blue.',
    });
    const [system, tools, user] = messagesOf(model, 0);
    expect(textOf(system as never)).toBe(guidance);
    const parts = (tools?.parts ?? []).flatMap((p) => (p.kind === 'tool-definition' ? [p] : []));
    expect(parts.filter((p) => !p.deferred).map((p) => p.definition.tool.id)).toEqual(
      toolkit.coreDefinitions().map((d) => d.tool.id),
    );
    expect(parts.filter((p) => p.deferred).map((p) => p.definition.tool.id)).toEqual(['archive', 'load_map']);
    expect(textOf(user as never)).toBe(
      'Memory: likes blue.\n\nCounter at 0.\n\nFocus: the display\n\nNotes since your last turn:\n- The person reset the counter.\n\nRequest: hello',
    );
  });

  it('a later turn sends only its request; the tools can all start loaded', async () => {
    const { session, model } = fixture([{ text: 'One.' }, { text: 'Two.' }]);
    await session.turn({ prompt: 'first' });
    await session.turn({ prompt: 'second' });
    const second = messagesOf(model, 1);
    expect(second.map((m) => m.role)).toEqual(['system', 'system', 'user', 'assistant', 'user']);
    expect(textOf(second.at(-1) as never)).toBe('Request: second');
    const all = fixture([{ text: 'Hi.' }], { startTools: 'all' });
    await all.session.turn({ prompt: 'hello' });
    const defs = (messagesOf(all.model, 0)[1]?.parts ?? []).flatMap((p) =>
      p.kind === 'tool-definition' ? [p.deferred] : [],
    );
    expect(defs.every((deferred) => !deferred)).toBe(true);
    expect(defs).toHaveLength(all.toolkit.allDefinitions().length);
  });

  it('a later session finds whether the guidance changed by its label, not by how its text ends', async () => {
    const storage = createMemoryStorage();
    const first = fixture([{ text: 'One.' }], { storage, harnessId: 'h', guidance: 'Be careful. Be brief.' });
    await first.session.turn({ prompt: 'first' });
    first.session.dispose();
    // The new guidance is how the old one ends: the text alone would call it unchanged.
    const second = fixture([{ text: 'Two.' }, { text: 'Three.' }], { storage, harnessId: 'h', guidance: 'Be brief.' });
    await second.session.turn({ prompt: 'second' });
    const sent = messagesOf(second.model, 0);
    const update = sent.find((m) => m.role === 'user' && textOf(m).startsWith(GUIDANCE_UPDATE));
    expect(textOf(update as never)).toBe(`${GUIDANCE_UPDATE}\n\nBe brief.`);
    // Once given, it is current: the next turn sends only its request.
    await second.session.turn({ prompt: 'third' });
    const third = messagesOf(second.model, 1);
    expect(third.filter((m) => textOf(m).startsWith(GUIDANCE_UPDATE))).toHaveLength(1);
  });

  it('forget() empties the conversation: the next turn is a first turn again', async () => {
    const { session, model } = fixture([{ text: 'One.' }, { text: 'Two.' }, { text: 'Three.' }]);
    await session.turn({ prompt: 'first' });
    expect(await session.isFresh()).toBe(false);
    await session.forget();
    expect(await session.isFresh()).toBe(true);
    expect(await headOf(session)).toEqual([]);
    await session.turn({ prompt: 'second' });
    expect(messagesOf(model, 1).map((m) => m.role)).toEqual(['system', 'system', 'user']);
    expect(resultsIn(await headOf(session))).toEqual([]);
    expect(textOf(messagesOf(model, 1).at(-1) as never)).toBe('Request: second');
  });
});
