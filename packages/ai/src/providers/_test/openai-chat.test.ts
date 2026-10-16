/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message, ModelVisiblePart } from '../../messages/index.js';
import { createOpenAIModelAdapter } from '../index.js';
import {
  answer,
  assertDigestMatches,
  collectFrames,
  lookupTool,
  openaiSse,
  payloadOf,
  project,
  recordingFetch,
  resultOf,
} from './support.js';

const secret = 'fixture-endpoint-secret';
const baseURL = 'https://llm.example.test/v1';
const limits = { 'kimi-k3': { contextWindowTokens: 262_144, maxOutputTokens: 16_384 } };
const create = (fetch?: typeof globalThis.fetch, api: 'chat' | 'responses' = 'chat') =>
  createOpenAIModelAdapter({
    resolveApiKey: () => secret,
    baseURL,
    api,
    provider: 'lab',
    label: 'Lab',
    modelLimits: limits,
    body: { cache: { 'no-cache': true } },
    ...(fetch === undefined ? {} : { fetch }),
  });
const input = (extras: { messages?: Message[]; parameters?: object } = {}) => ({
  configuration: {
    target: { provider: 'lab', model: 'kimi-k3' },
    parameters: { maxOutputTokens: 128, ...extras.parameters },
  },
  messages: extras.messages ?? [{ role: 'user' as const, parts: [{ kind: 'text' as const, text: 'Look up 1abc.' }] }],
  tools: [lookupTool],
});
const file = (mediaType: string): ModelVisiblePart => ({
  kind: 'reference',
  reference: { uri: `data:${mediaType};base64,${btoa('fixture')}`, mediaType },
});

/** A chat-completions stream: reasoning, then text or one tool call, then usage. */
function chatSse(options: { reasoning?: string; text?: string; call?: { name: string; input: string } }): string {
  const chunk = (choice: object) => ({
    id: 'chat_fixture',
    object: 'chat.completion.chunk',
    created: 1,
    model: 'kimi-k3',
    ...choice,
  });
  const delta = (d: object, finish: string | null = null) =>
    chunk({ choices: [{ index: 0, delta: d, finish_reason: finish }] });
  const events: object[] = [];
  if (options.reasoning) events.push(delta({ role: 'assistant', reasoning_content: options.reasoning }));
  if (options.text) events.push(delta({ content: options.text }));
  if (options.call) {
    const half = Math.floor(options.call.input.length / 2);
    events.push(
      delta({
        tool_calls: [
          {
            index: 0,
            id: 'call_fixture',
            type: 'function',
            function: { name: options.call.name, arguments: options.call.input.slice(0, half) },
          },
        ],
      }),
      delta({ tool_calls: [{ index: 0, function: { arguments: options.call.input.slice(half) } }] }),
    );
  }
  events.push(
    delta({}, options.call ? 'tool_calls' : 'stop'),
    chunk({
      choices: [],
      usage: {
        prompt_tokens: 100,
        completion_tokens: 20,
        prompt_tokens_details: { cached_tokens: 64 },
        completion_tokens_details: { reasoning_tokens: 5 },
      },
    }),
  );
  return `${events.map((e) => `data: ${JSON.stringify(e)}\n\n`).join('')}data: [DONE]\n\n`;
}

describe('OpenAI adapter, chat completions at another endpoint', () => {
  it('knows its models only by the limits given, has no tool additions or cache key, and checks the provider', async () => {
    const adapter = create();
    expect(await adapter.capabilitiesFor(input().configuration)).toMatchObject({
      contextWindowTokens: 262_144,
      appendsTools: false,
      cacheScope: false,
    });
    const luna = { target: { provider: 'lab', model: 'gpt-6-luna' }, parameters: {} };
    await expect(adapter.capabilitiesFor(luna)).rejects.toThrow(/explicit token limits/);
    const foreign = { target: { provider: 'openai', model: 'kimi-k3' }, parameters: {} };
    await expect(adapter.capabilitiesFor(foreign)).rejects.toThrow(/cannot project provider openai/);
  });

  it("streams reasoning and a call, sends the server's fields, and sends the reasoning back with the call", async () => {
    const recorded = recordingFetch(
      chatSse({ reasoning: 'look it up', call: { name: 'fixture_lookup', input: '{"id":"x"}' } }),
    );
    const adapter = create(recorded.fetch);
    const request = await project(adapter, input());
    await assertDigestMatches(request);
    const frames = await collectFrames(adapter, request);
    expect(frames).toContainEqual({ kind: 'reasoning-delta', text: 'look it up' });
    const result = resultOf(frames);
    expect(result.finish).toBe('tool-calls');
    expect(result.usage).toMatchObject({ inputTokens: 100, cacheReadInputTokens: 64, reasoningOutputTokens: 5 });
    expect(result.message.parts).toContainEqual({
      kind: 'tool-call',
      requestId: 'call_fixture',
      tool: lookupTool.tool,
      arguments: { id: 'x' },
      providerOptions: { chat: { text: 'look it up' } },
    });
    const first = recorded.requests[0];
    expect(first?.url).toBe(`${baseURL}/chat/completions`);
    expect(first?.headers.get('authorization')).toBe(`Bearer ${secret}`);
    expect(JSON.stringify(request.payload)).not.toContain(secret);
    const sent = JSON.parse(first?.body ?? '{}');
    expect(sent).toMatchObject({ model: 'kimi-k3', cache: { 'no-cache': true }, stream: true, max_tokens: 128 });
    expect(sent.tools.map((t: { function: { name: string } }) => t.function.name)).toEqual(['fixture_lookup']);
    expect(sent.reasoning_effort).toBeUndefined();

    const messages: Message[] = [...input().messages, result.message, answer('call_fixture')];
    await collectFrames(adapter, JSON.parse(JSON.stringify(await project(adapter, input({ messages })))));
    const replay = JSON.parse(recorded.requests[1]?.body ?? '{}').messages;
    expect(replay.slice(1)).toMatchObject([
      {
        role: 'assistant',
        reasoning_content: 'look it up',
        tool_calls: [{ id: 'call_fixture', function: { name: 'fixture_lookup', arguments: '{"id":"x"}' } }],
      },
      { role: 'tool', tool_call_id: 'call_fixture', content: 'found' },
    ]);
  });

  it('streams a text answer', async () => {
    const recorded = recordingFetch(chatSse({ text: 'PONG' }));
    const adapter = create(recorded.fetch);
    const result = resultOf(await collectFrames(adapter, await project(adapter, input())));
    expect(result).toMatchObject({ finish: 'complete', message: { parts: [{ kind: 'text', text: 'PONG' }] } });
  });

  it('sends an effort as reasoning_effort and refuses a token budget', async () => {
    const recorded = recordingFetch(chatSse({ text: 'ok' }));
    const adapter = create(recorded.fetch);
    const effort = { reasoning: { mode: 'adaptive', effort: 'low' } };
    await collectFrames(adapter, await project(adapter, input({ parameters: effort })));
    expect(JSON.parse(recorded.requests[0]?.body ?? '{}')).toMatchObject({
      reasoning_effort: 'low',
      cache: { 'no-cache': true },
    });
    const budget = { reasoning: { mode: 'budget', tokens: 1000 } };
    await expect(project(adapter, input({ parameters: budget }))).rejects.toThrow(/not a token budget/);
  });

  it('takes text files and images from the person, and refuses other attachments', async () => {
    const adapter = create();
    const user = (part: ModelVisiblePart): Message[] => [{ role: 'user', parts: [part] }];
    const image = payloadOf(await project(adapter, input({ messages: user(file('image/png')) })));
    expect(JSON.stringify(image.messages)).toContain('"mediaType":"image/png"');
    const text = payloadOf(await project(adapter, input({ messages: user(file('text/plain')) })));
    expect(JSON.stringify(text.messages)).toContain('"text":"fixture"');
    await expect(project(adapter, input({ messages: user(file('application/pdf')) }))).rejects.toThrow(
      /does not support inline application\/pdf/,
    );
    const call: Message = {
      role: 'assistant',
      parts: [{ kind: 'tool-call', requestId: 'c1', tool: lookupTool.tool, arguments: { id: 'x' } }],
    };
    const imageResult: Message = {
      role: 'user',
      parts: [
        {
          kind: 'tool-result',
          requestId: 'c1',
          tool: lookupTool.tool,
          status: 'succeeded',
          parts: [file('image/png')],
        },
      ],
    };
    await expect(
      project(adapter, input({ messages: [...user(file('text/plain')), call, imageResult] })),
    ).rejects.toThrow(/tool-result/);
  });

  it("speaks Responses at another endpoint when asked, without OpenAI's tool additions", async () => {
    const recorded = recordingFetch(openaiSse('PONG'));
    const adapter = create(recorded.fetch, 'responses');
    expect(await adapter.capabilitiesFor(input().configuration)).toMatchObject({
      appendsTools: false,
      cacheScope: true,
    });
    await collectFrames(adapter, await project(adapter, input()));
    expect(recorded.requests[0]?.url).toBe(`${baseURL}/responses`);
  });
});
