/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message, ModelConfiguration } from '@molstar/ai-agent';
import { createMemoryStorage } from '@molstar/ai-components';
import { createSession, type Session } from '@molstar/ai-conformance';
import { createToolBindings, createToolkit, NestedTools, RepeatedReads, toolFor } from '../index.js';
import { type Toy, toy } from './helpers.js';

/** Whole turns: a real harness and runner, the scripted model, the tools bound as a host binds them. */
const fixture = { name: 'memory', storage: () => createMemoryStorage() };
const none = { kind: 'object', fields: {} } as const;
const tool = toolFor<Toy, { name: string }>();

const loadedIn = (messages: Message[]) =>
  messages.flatMap((m) =>
    m.role === 'system'
      ? m.parts.flatMap((p) => (p.kind === 'tool-definition' && !p.deferred ? [p.definition.tool.id] : []))
      : [],
  );

async function setup(api: Toy, saved: string[] = []) {
  const get = tool({
    id: 'get_state',
    description: 'Reads the value.',
    input: none,
    modules: ['state'],
    read: true,
    run: ({ api }) => ({ ok: true, did: `Value ${api.state.value}.` }),
  });
  const set = tool({
    id: 'set_value',
    description: 'Bumps the value.',
    input: { kind: 'object', fields: { by: { kind: 'integer', optional: true } } },
    modules: ['state'],
    run: ({ args, api }) => {
      api.state.bump(args.by ?? 1);
      return { ok: true, did: 'Bumped.', changed: ['value'] };
    },
  });
  const load = tool({
    id: 'load_volume',
    description: 'Loads a map.',
    input: none,
    modules: ['log'],
    idempotent: true,
    run: ({ api, call }) => {
      api.log.lines.push(`load ${call.key} ${call.prior}`);
      return { ok: true, did: 'Loaded.' };
    },
  });
  const inner = tool({
    id: 'inner',
    description: 'The helper’s tool.',
    input: none,
    modules: ['state'],
    run: ({ api }) => {
      api.state.bump(10);
      return { ok: true, did: 'Inner done.' };
    },
  });
  const helper = tool({
    id: 'helper',
    description: 'Runs a helper.',
    input: none,
    modules: [],
    run: async ({ nest }) => {
      const r = await nest?.({
        messages: [{ role: 'user', parts: [{ kind: 'text', text: 'Help.' }] }],
        tools: [inner],
        api,
        model: { name: 'scripted-1' },
        limits: {},
      });
      return { ok: true, did: `Helper ${r?.status}.` };
    },
  });
  const kit = createToolkit<Toy>({
    tools: [get, set, load, helper],
    core: ['get_state', 'set_value', 'helper'],
    groups: [
      {
        id: 'volumes',
        title: 'Maps',
        summary: 'maps',
        guidance: 'Skill volumes: load a map.',
        tools: ['load_volume'],
      },
    ],
  });
  const configure = (choice: { name: string }): ModelConfiguration => ({
    target: { provider: 'scripted', model: choice.name },
    parameters: {},
  });
  const s: Session = await createSession(fixture, {
    model: [],
    tools: Object.fromEntries(kit.allDefinitions().map((definition) => [definition.tool.id, { definition }])),
    declare: kit.coreDefinitions().map((d) => d.tool.id),
    deferred: ['load_volume'],
  });
  api.disclosure = {
    loaded: async () => loadedIn(await s.main()),
    skills: async () => kit.groupsIn(await s.main()),
  };
  const bindings = createToolBindings(kit.tools, api, {
    repeats: new RepeatedReads({ version: () => api.state.version }),
    nested: new NestedTools<Toy, { name: string }>(configure),
    afterWrite: async ({ toolId }) => void saved.push(toolId),
  });
  return { s, bindings, kit };
}

const toolText = async (s: Session) => {
  const out: string[] = [];
  for (const m of await s.main())
    if (m.role === 'user')
      for (const p of m.parts)
        if (p.kind === 'tool-result') out.push(p.parts.map((x) => (x.kind === 'text' ? x.text : '')).join(''));
  return out;
};

describe('whole turns over the harness (tools-and-session.md §1, §3)', () => {
  it('a skill is disclosed, its tool then called, and the host’s write hook ran only for what changed', async () => {
    const api = toy();
    const saved: string[] = [];
    const { s, bindings } = await setup(api, saved);
    s.model.push(
      { calls: [{ tool: 'disclose_skills', arguments: { skills: ['volumes'] } }] },
      {
        calls: [
          { tool: 'load_volume', arguments: {} },
          { tool: 'set_value', arguments: { by: '2' } },
        ],
      },
      { calls: [{ tool: 'disclose_skills', arguments: { skills: ['volumes'] } }] },
      { text: 'Done.' },
    );
    const { exit } = await s.send('show the density', { runner: s.runner({ tools: bindings }) });
    expect(exit?.reason).toBe('released');
    expect((await s.turns()).map((t) => t.status)).toEqual(['completed']);
    expect(api.log.lines).toHaveLength(1);
    expect(api.log.lines[0]).toMatch(/^load \S+ none$/);
    expect(api.state.value).toBe(2);
    // disclose_skills changes nothing the host keeps: only the write ran the hook, and load_volume (a write too).
    expect(saved).toEqual(['load_volume', 'set_value']);
    const texts = await toolText(s);
    expect(texts[0]).toContain('Skill volumes loaded; its tools: load_volume.');
    expect(texts).toContain('Loaded.');
    expect(texts).toContain('Bumped.\nChanged: value.');
    expect(texts.at(-1)).toContain('Skill volumes is already loaded');
    expect(loadedIn(await s.main())).toContain('load_volume');
  });

  it('a read answered again in one line, and a refused call recorded as before any change', async () => {
    const api = toy();
    const { s, bindings } = await setup(api);
    s.model.push(
      { calls: [{ tool: 'get_state', arguments: {} }] },
      { calls: [{ tool: 'get_state', arguments: {} }] },
      { calls: [{ tool: 'set_value', arguments: { by: 'many' } }] },
      { text: 'Done.' },
    );
    await s.send('read twice', { runner: s.runner({ tools: bindings }) });
    const texts = await toolText(s);
    expect(texts[0]).toBe('Value 0.');
    expect(texts[1]).toContain('Unchanged: nothing changed since your earlier get_state call');
    expect(texts[2]).toBe('by: expected number (got "many").');
    // The failure is recorded as before any change (Stories never said, so every failure was `unknown`).
    const failures = (await s.events()).flatMap((e) =>
      e.kind === 'step' ? e.calls.flatMap((c) => c.attempts.flatMap((a) => (a.failure ? [a.failure.effect] : []))) : [],
    );
    expect(failures).toEqual(['none']);
    expect(api.state.value).toBe(0);
  });

  it('a tool runs a helper under its call, whose tool the helper’s model calls', async () => {
    const api = toy();
    const { s, bindings } = await setup(api);
    s.model.push(
      { calls: [{ tool: 'helper', arguments: {} }] },
      { calls: [{ tool: 'inner', arguments: {} }] },
      { text: 'Helper is done.' },
      { text: 'All done.' },
    );
    await s.send('use the helper', { runner: s.runner({ tools: bindings }) });
    expect((await s.turns()).map((t) => t.status)).toEqual(['completed']);
    expect(api.state.value).toBe(10);
    expect(await toolText(s)).toEqual(['Helper completed.']);
  });
});
