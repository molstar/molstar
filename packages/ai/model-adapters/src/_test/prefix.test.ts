/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { JsonObject, Message, ToolDefinition } from '@molstar/ai-agent';
import { UnrepresentableInput } from '@molstar/ai-agent';
import { createAnthropicModelAdapter, createGoogleModelAdapter, createOpenAIModelAdapter } from '../index.js';
import {
  answer,
  anthropicSse,
  collectFrames,
  declares,
  googleSse,
  lookupTool,
  payloadOf,
  project,
  projectInput,
  recordingFetch,
} from './support.js';

const guide: Message = { role: 'system', parts: [{ kind: 'text', text: 'Use the application-provided tools.' }] };
const user: Message = { role: 'user', parts: [{ kind: 'text', text: 'Find a value.' }] };
const tail: Message = { role: 'user', parts: [{ kind: 'text', text: 'Explain the result.' }] };
const other: ToolDefinition = {
  ...lookupTool,
  tool: { id: 'second.lookup' },
  name: 'second.lookup',
  description: 'Second tool with a complete schema.',
};
const sonnet = projectInput('anthropic', 'claude-sonnet-5-5');

function body(recorded: ReturnType<typeof recordingFetch>, index: number): JsonObject {
  return JSON.parse(recorded.requests[index]?.body ?? '{}');
}
function content(message: JsonObject | undefined): JsonObject[] {
  return Array.isArray(message?.content) ? (message.content as JsonObject[]) : [];
}
function withoutCache(value: unknown): unknown {
  if (Array.isArray(value)) return value.map(withoutCache);
  if (value === null || typeof value !== 'object') return value;
  return Object.fromEntries(
    Object.entries(value)
      .filter(([key]) => key !== 'cache_control')
      .map(([key, item]) => [key, withoutCache(item)]),
  );
}

describe('provider prompt prefix mechanics', () => {
  // harness.md §11: the request's parts (tools, system, each message), so a trace can say where a prefix broke.
  it.each(['anthropic', 'google', 'openai'] as const)(
    'gives %s requests their parts: an appended message extends them, other tools change the first',
    async (provider) => {
      const recorded = recordingFetch(anthropicSse('ok'));
      const create = {
        anthropic: createAnthropicModelAdapter,
        google: createGoogleModelAdapter,
        openai: createOpenAIModelAdapter,
      };
      const adapter = create[provider]({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
      const model = { anthropic: 'claude-haiku-4-5', google: 'gemini-3.8-flash', openai: 'gpt-6-luna' }[provider];
      const at = (messages: Message[], tools: ToolDefinition[]) =>
        project(adapter, projectInput(provider, model, { messages, tools }));
      const first = (await at([guide, user], [lookupTool])).segments ?? [];
      const appended = (await at([guide, user, tail], [lookupTool])).segments ?? [];
      const retooled = (await at([guide, user], [lookupTool, other])).segments ?? [];
      expect(first.length).toBeGreaterThan(2);
      expect(appended.slice(0, first.length)).toEqual(first);
      expect(appended.length).toBeGreaterThan(first.length);
      // Other tools break the prefix at the tools (OpenAI also loads a tool given late as a message).
      expect(retooled[0]).not.toBe(first[0]);
      if (provider !== 'openai') expect(retooled.slice(1)).toEqual(first.slice(1));
    },
  );

  it.each(['anthropic', 'google'] as const)(
    'preserves %s prompt content when a user message is appended',
    async (provider) => {
      const recorded = recordingFetch(provider === 'anthropic' ? anthropicSse('ok') : googleSse('ok'));
      const adapter =
        provider === 'anthropic'
          ? createAnthropicModelAdapter({ resolveApiKey: () => 'fixture', fetch: recorded.fetch })
          : createGoogleModelAdapter({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
      const model = provider === 'anthropic' ? 'claude-haiku-4-5' : 'gemini-3.8-flash';
      const request = await project(adapter, projectInput(provider, model, { messages: [guide, user] }));
      const appended = await project(adapter, projectInput(provider, model, { messages: [guide, user, tail] }));
      await collectFrames(adapter, request);
      await collectFrames(adapter, appended);
      const first = body(recorded, 0);
      const second = body(recorded, 1);
      if (provider === 'anthropic') {
        expect(withoutCache(second.system)).toEqual(withoutCache(first.system));
        expect(withoutCache(second.tools)).toEqual(withoutCache(first.tools));
        // Anthropic merges adjacent user turns; the old content stays an exact prefix inside the merged turn.
        expect((second.messages as JsonObject[])[0]?.content).toEqual([
          ...content((first.messages as JsonObject[])[0]),
          { type: 'text', text: 'Explain the result.' },
        ]);
      } else {
        expect(second.systemInstruction).toEqual(first.systemInstruction);
        expect(second.tools).toEqual(first.tools);
        expect((second.contents as unknown[]).slice(0, (first.contents as unknown[]).length)).toEqual(first.contents);
      }
    },
  );

  it.each([
    ['anthropic', createAnthropicModelAdapter, 'claude-sonnet-5-5'],
    ['google', createGoogleModelAdapter, 'gemini-3.8-flash'],
    ['openai', createOpenAIModelAdapter, 'gpt-6-luna'],
  ] as const)(
    'keeps the %s payload prefix byte-stable when only messages are appended',
    async (provider, create, model) => {
      const adapter = create({ resolveApiKey: () => 'fixture' });
      const base = [declares(lookupTool), guide, user];
      const next = [...base, { role: 'assistant' as const, parts: [{ kind: 'text' as const, text: 'ok' }] }, tail];
      const input = (messages: Message[]) => projectInput(provider, model, { messages, tools: [lookupTool] });
      const first = payloadOf(await project(adapter, input(base)));
      const second = payloadOf(await project(adapter, input(next)));
      expect(second.system).toBe(first.system);
      expect(second.tools).toEqual(first.tools);
      expect((second.messages as unknown[]).slice(0, (first.messages as unknown[]).length)).toEqual(first.messages);
    },
  );

  it('appends native full-schema disclosure without changing the cache-relevant earlier prompt', async () => {
    const recorded = recordingFetch(anthropicSse('ok'));
    const adapter = createAnthropicModelAdapter({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
    const head: Message[] = [
      { role: 'system', parts: [...guide.parts, { kind: 'tool-definition', definition: lookupTool, deferred: false }] },
      declares(other, true),
      user,
    ];
    const first = await project(adapter, { ...sonnet, messages: head, tools: [lookupTool] });
    const second = await project(adapter, {
      ...sonnet,
      messages: [...head, declares(other)],
      tools: [lookupTool, other],
    });
    await collectFrames(adapter, first);
    await collectFrames(adapter, second);
    const initial = body(recorded, 0);
    const extended = body(recorded, 1);
    expect(initial.tools).toEqual(extended.tools);
    expect(initial.system).toEqual(extended.system);
    expect((extended.messages as unknown[]).slice(0, (initial.messages as unknown[]).length)).toEqual(initial.messages);
    expect((extended.messages as JsonObject[]).at(-1)).toEqual({
      role: 'system',
      content: [{ type: 'tool_addition', tool: { type: 'tool_reference', name: 'second_lookup' } }],
    });
    const catalog = extended.tools as JsonObject[];
    expect(catalog[1]?.defer_loading).toBe(true);
    expect(catalog[1]?.input_schema).toEqual(other.inputSchema);
    expect(recorded.requests[1]?.headers.get('anthropic-beta')).toContain('mid-conversation-tool-changes-2026-07-01');
    expect(second.inputTokens).toBeGreaterThan(first.inputTokens);
  });

  it('supports legal in-place system guidance and rejects unrepresentable placements', async () => {
    const recorded = recordingFetch(anthropicSse('ok'));
    const adapter = createAnthropicModelAdapter({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
    const update: Message = { role: 'system', parts: [{ kind: 'text', text: 'New operator constraint.' }] };
    const input = { ...sonnet, tools: [] };
    await collectFrames(adapter, await project(adapter, { ...input, messages: [guide, user, update] }));
    expect((body(recorded, 0).messages as JsonObject[]).at(-1)).toEqual({
      role: 'system',
      content: [{ type: 'text', text: 'New operator constraint.' }],
    });
    const boundary = await project(adapter, { ...input, messages: [guide, user, update, user] }).catch(
      (e: unknown) => e,
    );
    expect(boundary).toBeInstanceOf(UnrepresentableInput);
    expect((boundary as Error).message).toContain('boundary');
    const noInitial = await project(adapter, { ...input, messages: [user, update] }).catch((e: unknown) => e);
    expect(noInitial).toBeInstanceOf(UnrepresentableInput);
    expect((noInitial as Error).message).toContain('initial system prompt');
  });

  it('keeps signed thinking and attachment results adjacent before an ordered load', async () => {
    const recorded = recordingFetch(anthropicSse('ok'));
    const adapter = createAnthropicModelAdapter({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
    const response: Message = {
      role: 'assistant',
      parts: [
        {
          kind: 'tool-call',
          requestId: 'call',
          tool: lookupTool.tool,
          arguments: { id: 'value' },
          providerOptions: { anthropic: { thinking: [{ text: 'reason', signature: 'signed' }] } },
        },
      ],
    };
    const result: Message = {
      role: 'user',
      parts: [
        {
          kind: 'tool-result',
          requestId: 'call',
          tool: lookupTool.tool,
          status: 'succeeded',
          parts: [
            { kind: 'text', text: 'result' },
            { kind: 'reference', reference: { uri: 'data:image/png;base64,YQ==' } },
          ],
        },
      ],
    };
    await collectFrames(
      adapter,
      await project(adapter, {
        ...sonnet,
        messages: [guide, declares(lookupTool), user, response, result, declares(other)],
        tools: [lookupTool, other],
      }),
    );
    const messages = body(recorded, 0).messages as JsonObject[];
    expect(content(messages[1])[0]).toEqual({ type: 'thinking', thinking: 'reason', signature: 'signed' });
    expect(messages.map((message) => message.role)).toEqual(['user', 'assistant', 'user', 'system']);
    expect(content(content(messages[2])[0]).map((part) => part.type)).toEqual(['text', 'image']);
    expect(content(messages[3])[0]?.type).toBe('tool_addition');
  });

  it('rejects version histories that native reference loading cannot reproduce', async () => {
    const adapter = createAnthropicModelAdapter({ resolveApiKey: () => 'unused' });
    const updated = {
      ...lookupTool,
      description: 'New schema version',
      inputSchema: { type: 'object', properties: { id: { enum: ['new'] } } },
    };
    await expect(
      project(adapter, {
        ...sonnet,
        messages: [declares(lookupTool), guide, user, declares(updated)],
        tools: [updated],
      }),
    ).rejects.toThrow('use collected projection or a new context head');
  });

  it('keeps cache activation adapter-owned, bounded, and optional', async () => {
    for (const cache of [true, false]) {
      const recorded = recordingFetch(anthropicSse('ok'));
      const adapter = createAnthropicModelAdapter({ resolveApiKey: () => 'fixture', fetch: recorded.fetch, cache });
      await collectFrames(adapter, await project(adapter, { ...sonnet, messages: [guide, user] }));
      const wire = body(recorded, 0);
      expect(JSON.stringify(wire).match(/"cache_control"/g)?.length ?? 0).toBe(cache ? 3 : 0);
      if (cache) expect(wire.cache_control).toEqual({ type: 'ephemeral', ttl: '5m' });
    }
  });

  it('negotiates target-specific tool appending and rejects capabilities a target cannot honor', async () => {
    const claude = createAnthropicModelAdapter({ resolveApiKey: () => 'unused' });
    expect((await claude.capabilitiesFor(sonnet.configuration)).appendsTools).toBe(true);
    const old = projectInput('anthropic', 'claude-haiku-4-5', { messages: [user], tools: [] });
    const oldCapabilities = await claude.capabilitiesFor(old.configuration);
    expect(oldCapabilities.appendsTools).toBe(false);
    await expect(project(claude, old, { ...oldCapabilities, appendsTools: true })).rejects.toThrow('does not support');
    const google = createGoogleModelAdapter({ resolveApiKey: () => 'unused' });
    const gemini = projectInput('google', 'gemini-3.8-flash');
    const geminiCapabilities = await google.capabilitiesFor(gemini.configuration);
    await expect(project(google, gemini, { ...geminiCapabilities, appendsTools: true })).rejects.toThrow(
      'does not support',
    );
  });

  it('answers a call and loads a tool in the same history, in order', async () => {
    const adapter = createAnthropicModelAdapter({ resolveApiKey: () => 'unused' });
    const call: Message = {
      role: 'assistant',
      parts: [{ kind: 'tool-call', requestId: 'c1', tool: lookupTool.tool, arguments: {} }],
    };
    const request = await project(adapter, {
      ...sonnet,
      messages: [guide, declares(lookupTool), user, call, answer('c1'), declares(other)],
      tools: [lookupTool, other],
    });
    expect((payloadOf(request).messages as JsonObject[]).map((message) => message.role)).toEqual([
      'user',
      'assistant',
      'tool',
      'system',
    ]);
  });
});
