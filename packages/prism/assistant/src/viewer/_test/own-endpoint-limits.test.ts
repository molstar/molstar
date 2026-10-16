/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AgentSession, createToolkit } from '@molstar/ai';
import { describe, expect, it } from 'vitest';
import type { KeyScope, KeyStore } from '../../keys.js';
import { checkLimits, OwnLimits } from '../../model.js';
import { createAISettings, type SettingsPrefs } from '../settings-model.js';

// An own OpenAI-compatible endpoint runs turns over the real provider (agent.md §4): its models' limits are the
// person's, else 32,768 tokens of context and 4,096 of output. The endpoint is a recorded fetch: no network.

function memoryKeys(): KeyStore {
  const stored = new Map<string, string>();
  const scopes = new Map<string, KeyScope>();
  return {
    get: async (p) => stored.get(p),
    async set(p, key, scope) {
      stored.set(p, key);
      scopes.set(p, scope);
    },
    async remove(p) {
      stored.delete(p);
      scopes.delete(p);
    },
    list: async () => [...scopes].map(([provider, scope]) => ({ provider, scope })),
    protection: async () => 'per-browser',
  };
}

function memoryPrefs(): SettingsPrefs {
  const values = new Map<string, string>();
  return { get: (k) => values.get(k), set: (k, v) => void values.set(k, v) };
}

/** A chat-completions stream that says `text`, as an OpenAI-compatible server sends it. */
function chatSse(text: string): string {
  const chunk = (fields: object) => ({
    id: 'chat_fixture',
    object: 'chat.completion.chunk',
    created: 1,
    model: 'llama',
    ...fields,
  });
  const events = [
    chunk({ choices: [{ index: 0, delta: { role: 'assistant', content: text }, finish_reason: null }] }),
    chunk({ choices: [{ index: 0, delta: {}, finish_reason: 'stop' }] }),
    chunk({ choices: [], usage: { prompt_tokens: 100, completion_tokens: 5 } }),
  ];
  return `${events.map((e) => `data: ${JSON.stringify(e)}\n\n`).join('')}data: [DONE]\n\n`;
}

/** The endpoint: every request recorded, every answer `text`. */
function endpoint(text: string) {
  const requests: { url: string; body: { max_tokens?: number; model?: string } }[] = [];
  const fetch: typeof globalThis.fetch = async (input, init) => {
    requests.push({ url: String(input), body: JSON.parse(typeof init?.body === 'string' ? init.body : '{}') });
    return new Response(chatSse(text), { status: 200, headers: { 'content-type': 'text/event-stream' } });
  };
  return { fetch, requests };
}

async function ownTurn(limits: { contextWindow?: number; maxOutputTokens?: number }) {
  const server = endpoint('Hello from the lab.');
  const settings = createAISettings(memoryKeys(), memoryPrefs(), { fetch: server.fetch });
  await settings.load();
  settings.addProvider({
    label: 'Lab',
    baseURL: 'http://localhost:11434/v1',
    api: 'chat',
    key: '',
    scope: 'browser',
    models: ['llama'],
    ...limits,
  });
  const setup = settings.setup.value;
  expect(setup).toEqual({ provider: 'own-lab', model: 'llama' });
  const session = new AgentSession({
    models: { providerFor: (id) => settings.providerFor(id), default: setup! },
    toolkit: createToolkit<Record<string, never>>({ tools: [], core: [], groups: [] }),
    api: {},
    guidance: 'You help.',
  });
  const result = await session.turn({ prompt: 'Hi.' });
  return { result, server, session };
}

describe('an own endpoint’s limits', () => {
  it('runs a turn with none given: 32,768 tokens of context, 4,096 of output', async () => {
    const { result, server, session } = await ownTurn({});
    expect(result.error).toBeUndefined();
    expect(result).toMatchObject({ status: 'completed', text: 'Hello from the lab.' });
    expect(server.requests.map((r) => [r.url, r.body.model, r.body.max_tokens])).toEqual([
      ['http://localhost:11434/v1/chat/completions', 'llama', OwnLimits.maxOutputTokens],
    ]);
    expect(session.contextWindow()).toBe(OwnLimits.contextWindowTokens);
  });

  it('runs a turn with the limits the person gave', async () => {
    const { result, server, session } = await ownTurn({ contextWindow: 131_072, maxOutputTokens: 8192 });
    expect(result).toMatchObject({ status: 'completed', text: 'Hello from the lab.' });
    expect(server.requests.map((r) => r.body.max_tokens)).toEqual([8192]);
    expect(session.contextWindow()).toBe(131_072);
  });

  it('refuses limits that are not positive whole numbers, or an output that fills the window', () => {
    expect(checkLimits({})).toBeUndefined();
    expect(checkLimits({ contextWindow: 65_536, maxOutputTokens: 8192 })).toBeUndefined();
    expect(checkLimits({ contextWindow: 0 })).toContain('whole numbers');
    expect(checkLimits({ maxOutputTokens: 1.5 })).toContain('whole numbers');
    expect(checkLimits({ maxOutputTokens: Number.NaN })).toContain('whole numbers');
    expect(checkLimits({ contextWindow: 4096 })).toContain('less than the context window');
    const reported: string[] = [];
    const settings = createAISettings(memoryKeys(), memoryPrefs(), { report: (t) => reported.push(t) });
    settings.addProvider({
      label: 'Lab',
      baseURL: 'http://localhost:11434/v1',
      api: 'chat',
      key: '',
      scope: 'browser',
      models: ['llama'],
      contextWindow: -1,
    });
    expect(reported).toEqual([
      'Lab was not added: The context window and the max output are whole numbers of tokens, above zero.',
    ]);
    expect(settings.providers.value.some((p) => p.id === 'own-lab')).toBe(false);
  });

  it('keeps the limits across a reload, and a key given later', async () => {
    const prefs = memoryPrefs();
    const keys = memoryKeys();
    const settings = createAISettings(keys, prefs);
    await settings.load();
    settings.addProvider({
      label: 'Lab',
      baseURL: 'http://localhost:11434/v1',
      api: 'chat',
      key: '',
      scope: 'browser',
      models: ['llama'],
      contextWindow: 65_536,
    });
    settings.setKey('own-lab', 'sk-lab', 'browser');
    expect(JSON.parse(prefs.get('providers') ?? '[]')).toEqual([
      expect.objectContaining({ id: 'own-lab', contextWindow: 65_536 }),
    ]);
  });
});
