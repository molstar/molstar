/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message, ToolDefinition } from '../../messages/index.js';
import { createAnthropicProvider, createGoogleProvider, createOpenAIProvider } from '../index.js';
import { declares, lookupTool, payloadOf, project, projectInput } from './support.js';

const guide: Message = { role: 'system', parts: [{ kind: 'text', text: 'Use the application-provided tools.' }] };
const user: Message = { role: 'user', parts: [{ kind: 'text', text: 'Find a value.' }] };
const reply: Message = { role: 'assistant', parts: [{ kind: 'text', text: 'ok' }] };
const tail: Message = { role: 'user', parts: [{ kind: 'text', text: 'Explain the result.' }] };
const other: ToolDefinition = {
  ...lookupTool,
  tool: { id: 'second.lookup' },
  name: 'second.lookup',
  description: 'Second tool with a complete schema.',
};

// Providers cache by prefix: an appended turn's request starts with the previous one, and its parts (tools, system,
// each message: harness.md §11) say so, so the trace can tell where a prefix broke.
describe('an appended turn', () => {
  it.each([
    ['anthropic', createAnthropicProvider, 'claude-sonnet-5-5'],
    ['anthropic', createAnthropicProvider, 'claude-haiku-4-5'],
    ['google', createGoogleProvider, 'gemini-3.8-flash'],
    ['openai', createOpenAIProvider, 'gpt-6-luna'],
  ] as const)('starts with the previous request: %s %s', async (provider, create, model) => {
    const adapter = create({ resolveApiKey: () => 'fixture' });
    const at = (messages: Message[], tools: ToolDefinition[] = [lookupTool]) =>
      project(adapter, projectInput(provider, model, { messages, tools }));
    const base = [declares(lookupTool), guide, user];
    const first = await at(base);
    const second = await at([...base, reply, tail]);
    const [a, b] = [payloadOf(first), payloadOf(second)];
    expect(b.system).toEqual(a.system);
    expect(b.tools).toEqual(a.tools);
    expect((b.messages as unknown[]).slice(0, (a.messages as unknown[]).length)).toEqual(a.messages);
    expect((first.segments ?? []).length).toBeGreaterThan(2);
    expect(second.segments?.slice(0, first.segments?.length)).toEqual(first.segments);
    // Other tools break the prefix at the tools (the appended layouts load a tool given late as a message).
    const retooled = (await at([guide, user], [lookupTool, other])).segments ?? [];
    const plain = (await at([guide, user])).segments ?? [];
    expect(retooled[0]).not.toBe(plain[0]);
    if (model === 'claude-haiku-4-5' || provider === 'google') expect(retooled.slice(1)).toEqual(plain.slice(1));
  });
});
