/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message, ModelAdapter } from '../../messages/index.js';
import { createAnthropicModelAdapter, createGoogleModelAdapter, createOpenAIModelAdapter } from '../index.js';
import { answer, lookupTool, project } from './support.js';

const key = { resolveApiKey: () => 'fixture' };
const adapters: readonly { name: string; provider: string; model: string; adapter: ModelAdapter }[] = [
  {
    name: 'anthropic native',
    provider: 'anthropic',
    model: 'claude-sonnet-5-5',
    adapter: createAnthropicModelAdapter(key),
  },
  {
    name: 'anthropic collected',
    provider: 'anthropic',
    model: 'claude-haiku-4-5',
    adapter: createAnthropicModelAdapter(key),
  },
  { name: 'openai responses', provider: 'openai', model: 'gpt-6-luna', adapter: createOpenAIModelAdapter(key) },
  {
    name: 'openai chat',
    provider: 'lab',
    model: 'lab-model',
    adapter: createOpenAIModelAdapter({
      ...key,
      api: 'chat',
      baseURL: 'https://lab.example/v1',
      provider: 'lab',
      modelLimits: { 'lab-model': { contextWindowTokens: 100_000, maxOutputTokens: 4_096 } },
    }),
  },
  { name: 'google', provider: 'google', model: 'gemini-3.8-flash', adapter: createGoogleModelAdapter(key) },
];

/** A conversation with every kind of entry a session labels. */
const labelled: Message[] = [
  { role: 'system', parts: [{ kind: 'text', text: 'Be brief.' }], label: { kind: 'session:guidance', key: 'g1' } },
  {
    role: 'system',
    parts: [{ kind: 'tool-definition', definition: lookupTool, deferred: false }],
    label: { kind: 'tools:definitions' },
  },
  {
    role: 'user',
    parts: [{ kind: 'text', text: 'State: 1abc.\n\nRequest: Look it up.' }],
    label: { kind: 'session:request', data: { prompt: 'Look it up.' } },
  },
  {
    role: 'assistant',
    parts: [{ kind: 'tool-call', requestId: 'r1', tool: lookupTool.tool, arguments: { id: '1abc' } }],
  },
  answer('r1'),
  { role: 'assistant', parts: [{ kind: 'text', text: 'Found.' }] },
  {
    role: 'user',
    parts: [{ kind: 'text', text: 'Request: And again.' }],
    label: { kind: 'session:request', data: { prompt: 'And again.' } },
  },
];
const unlabelled = labelled.map(({ role, parts }) => ({ role, parts }) as Message);

describe('labels are never sent', () => {
  it.each(adapters)('$name: a labelled conversation projects exactly as an unlabelled one', async (c) => {
    const input = (messages: Message[]) => ({
      configuration: { target: { provider: c.provider, model: c.model }, parameters: { maxOutputTokens: 128 } },
      messages,
      tools: [lookupTool],
    });
    const a = await project(c.adapter, input(labelled));
    const b = await project(c.adapter, input(unlabelled));
    expect(a.payload).toEqual(b.payload);
    expect(a.digest).toBe(b.digest);
    expect(a.inputTokens).toBe(b.inputTokens);
    expect(a.segments).toEqual(b.segments);
    expect(JSON.stringify(a.payload)).not.toContain('session:request');
  });
});
