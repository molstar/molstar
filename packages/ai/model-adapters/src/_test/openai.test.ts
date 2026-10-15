/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { digest, type Message, type ModelParameters, type ModelVisiblePart } from '@molstar/ai-agent';
import { createOpenAIModelAdapter } from '../index.js';
import {
  answer,
  assertDigestMatches,
  collectFrames,
  declares,
  lookupTool,
  openaiSse,
  payloadOf,
  project,
  projectInput,
  recordingFetch,
  resultOf,
} from './support.js';

const secret = 'fixture-openai-secret';
const input = (extras: Parameters<typeof projectInput>[2] = {}) => projectInput('openai', 'gpt-6-luna', extras);
const create = (fetch?: typeof globalThis.fetch) =>
  createOpenAIModelAdapter({ resolveApiKey: () => secret, ...(fetch === undefined ? {} : { fetch }) });
const file = (mediaType: string, text = 'fixture'): ModelVisiblePart => ({
  kind: 'reference',
  reference: { uri: `data:${mediaType};base64,${btoa(text)}`, mediaType },
});

describe('OpenAI Responses adapter', () => {
  it('resolves Luna limits and features, rejects unknown models, and validates override and output caps', async () => {
    const adapter = create();
    expect(await adapter.capabilitiesFor(input().configuration)).toMatchObject({
      contextWindowTokens: 1_050_000,
      maxOutputTokens: 128_000,
      appendsTools: true,
      cacheScope: true,
      outputLimitIncludesReasoning: true,
    });
    const unknown = { target: { provider: 'openai', model: 'unknown' }, parameters: {} };
    await expect(adapter.capabilitiesFor(unknown)).rejects.toThrow(/explicit token limits/);
    const custom = createOpenAIModelAdapter({
      resolveApiKey: () => secret,
      modelLimits: { gateway: { contextWindowTokens: 10_000, maxOutputTokens: 1024, defaultMaxOutputTokens: 256 } },
    });
    const config = { target: { provider: 'openai', model: 'gateway' }, parameters: {} };
    const gateway = await custom.capabilitiesFor(config);
    expect(gateway).toMatchObject({ appendsTools: false, cacheScope: true });
    expect(payloadOf(await project(custom, { configuration: config, messages: [], tools: [] })).maxOutputTokens).toBe(
      256,
    );
    await expect(
      project(custom, { configuration: { ...config, parameters: { maxOutputTokens: 1025 } }, messages: [], tools: [] }),
    ).rejects.toThrow(/model maximum/);
  });

  it('streams text and usage through stateless Responses without credential evidence', async () => {
    const recorded = recordingFetch(openaiSse('PONG'));
    const adapter = create(recorded.fetch);
    const request = await project(adapter, input({ tools: [] }));
    const frames = await collectFrames(adapter, request);
    await assertDigestMatches(request);
    expect(frames).toContainEqual({ kind: 'text-delta', text: 'PONG' });
    expect(resultOf(frames)).toMatchObject({
      finish: 'complete',
      usage: { inputTokens: 12, outputTokens: 11, cacheReadInputTokens: 4, reasoningOutputTokens: 3 },
    });
    const body = JSON.parse(recorded.requests[0]?.body ?? '{}');
    expect(recorded.requests[0]?.url).toMatch(/\/responses$/);
    expect(body).toMatchObject({ model: 'gpt-6-luna', store: false, max_output_tokens: 128, stream: true });
    expect(JSON.stringify(request)).not.toContain(secret);
    expect(JSON.stringify(frames)).not.toContain(secret);
    expect(recorded.requests[0]?.body).not.toContain(secret);
  });

  it.each(['none', 'low', 'medium', 'high', 'xhigh', 'max'])('maps reasoning effort %s', async (effort) => {
    const recorded = recordingFetch(openaiSse('done'));
    const adapter = create(recorded.fetch);
    const request = await project(
      adapter,
      input({
        configuration: {
          ...input().configuration,
          parameters: { maxOutputTokens: 128, reasoning: { mode: 'adaptive', effort } },
        },
      }),
    );
    await collectFrames(adapter, request);
    expect(JSON.parse(recorded.requests[0]?.body ?? '{}').reasoning.effort).toBe(effort);
  });

  it("maps explicit off to none instead of allowing the model's default reasoning", async () => {
    const recorded = recordingFetch(openaiSse('done'));
    const adapter = create(recorded.fetch);
    const request = await project(
      adapter,
      input({
        configuration: { ...input().configuration, parameters: { maxOutputTokens: 128, reasoning: { mode: 'off' } } },
      }),
    );
    await collectFrames(adapter, request);
    expect(JSON.parse(recorded.requests[0]?.body ?? '{}').reasoning.effort).toBe('none');
  });

  it.each([
    { reasoning: { mode: 'budget', tokens: 1024 } },
    { reasoning: { mode: 'adaptive', effort: 'minimal' } },
    { temperature: 0.5 },
    { topP: 0.5 },
  ] satisfies readonly ModelParameters[])(
    'rejects unsupported settings as a plain error before fetch %#',
    async (parameters) => {
      const recorded = recordingFetch(openaiSse('done'));
      const adapter = create(recorded.fetch);
      const failure = await project(adapter, input({ configuration: { ...input().configuration, parameters } })).catch(
        (error: unknown) => error,
      );
      expect((failure as Error).message).toMatch(/unsupported|not supported/);
      expect((failure as Error).name).toBe('Error');
      expect(recorded.requests).toHaveLength(0);
    },
  );

  it('retains encrypted reasoning and exact function-call identity for a subsequent request', async () => {
    const recorded = recordingFetch(
      openaiSse(
        '',
        { name: 'fixture_lookup', input: '{"id":"x"}' },
        { text: 'look up', encryptedContent: 'encrypted-fixture' },
      ),
    );
    const adapter = create(recorded.fetch);
    const result = resultOf(await collectFrames(adapter, await project(adapter, input())));
    expect(result.finish).toBe('tool-calls');
    expect(result.message.parts.find((part) => part.kind === 'tool-call')).toMatchObject({
      kind: 'tool-call',
      requestId: 'call_fixture',
      tool: lookupTool.tool,
      arguments: { id: 'x' },
      providerOptions: {
        openai: {
          itemId: 'fc_fixture',
          reasoning: [{ itemId: 'rs_fixture', text: 'look up', encryptedContent: 'encrypted-fixture' }],
        },
      },
    });
    const messages: Message[] = [...input().messages, result.message, answer('call_fixture')];
    const replay = await project(adapter, input({ messages }));
    await collectFrames(adapter, JSON.parse(JSON.stringify(replay)));
    const body = JSON.parse(recorded.requests[1]?.body ?? '{}');
    const reasonIndex = body.input.findIndex((item: { type: string }) => item.type === 'reasoning');
    expect(body.input[reasonIndex]).toMatchObject({
      id: 'rs_fixture',
      encrypted_content: 'encrypted-fixture',
      summary: [{ type: 'summary_text', text: 'look up' }],
    });
    expect(body.input[reasonIndex + 1]).toMatchObject({
      type: 'function_call',
      call_id: 'call_fixture',
      name: 'fixture_lookup',
      arguments: '{"id":"x"}',
    });
    expect(body.input[reasonIndex + 2]).toMatchObject({
      type: 'function_call_output',
      call_id: 'call_fixture',
      output: 'found',
    });
  });

  it('retains encrypted reasoning even when the provider supplies no visible summary', async () => {
    const wire = openaiSse(
      '',
      { name: 'fixture_lookup', input: '{"id":"x"}' },
      { text: '', encryptedContent: 'opaque-fixture' },
    )
      .split('\n')
      .filter((line) => !line.includes('response.reasoning_summary'))
      .join('\n');
    const adapter = create(recordingFetch(wire).fetch);
    expect(resultOf(await collectFrames(adapter, await project(adapter, input())))).toMatchObject({
      message: {
        parts: [
          {
            kind: 'tool-call',
            providerOptions: {
              openai: { reasoning: [{ itemId: 'rs_fixture', text: '', encryptedContent: 'opaque-fixture' }] },
            },
          },
        ],
      },
    });
  });

  it('keeps the SDK prompt prefix when a user message is appended', async () => {
    const recorded = recordingFetch(openaiSse('done'));
    const adapter = create(recorded.fetch);
    const guide: Message = { role: 'system', parts: [{ kind: 'text', text: 'policy' }] };
    const first = await project(adapter, input({ messages: [guide, ...input().messages] }));
    const tail: Message = { role: 'user', parts: [{ kind: 'text', text: 'tail' }] };
    const second = await project(adapter, input({ messages: [guide, ...input().messages, tail] }));
    await collectFrames(adapter, first);
    await collectFrames(adapter, second);
    const before = JSON.parse(recorded.requests[0]?.body ?? '{}');
    const after = JSON.parse(recorded.requests[1]?.body ?? '{}');
    expect(after.input.slice(0, before.input.length)).toEqual(before.input);
    expect(after.tools).toEqual(before.tools);
  });

  it('projects ordered user and tool-result images/documents and rejects unsupported media placements', async () => {
    const recorded = recordingFetch(openaiSse('done'));
    const adapter = create(recorded.fetch);
    const messages: Message[] = [
      {
        role: 'user',
        parts: [
          { kind: 'text', text: 'before' },
          file('image/png'),
          file('application/pdf'),
          file('text/markdown', '# hello'),
          file('application/vnd.openxmlformats-officedocument.wordprocessingml.document'),
        ],
      },
      {
        role: 'assistant',
        parts: [{ kind: 'tool-call', requestId: 'call-fixture', tool: lookupTool.tool, arguments: { id: 'x' } }],
      },
      {
        role: 'user',
        parts: [
          { kind: 'text', text: 'after' },
          file('image/jpeg'),
          {
            kind: 'tool-result',
            requestId: 'call-fixture',
            tool: lookupTool.tool,
            status: 'succeeded',
            parts: [file('image/webp'), file('application/pdf'), { kind: 'text', text: 'found' }],
          },
        ],
      },
    ];
    await collectFrames(adapter, await project(adapter, input({ messages, tools: [] })));
    const body = JSON.parse(recorded.requests[0]?.body ?? '{}');
    const types = (parts: { type: string }[]) => parts.map((part) => part.type);
    expect(types(body.input[0].content)).toEqual([
      'input_text',
      'input_image',
      'input_file',
      'input_text',
      'input_file',
    ]);
    expect(body.input[0].content[2].file_data).toBe(`data:application/pdf;base64,${btoa('fixture')}`);
    expect(body.input[2].type).toBe('function_call_output');
    expect(types(body.input[2].output)).toEqual(['input_image', 'input_file', 'input_text']);
    expect(types(body.input[3].content)).toEqual(['input_text', 'input_image']);
    for (const [role, part] of [
      ['user', file('audio/wav')],
      ['user', file('video/mp4')],
      ['user', file('application/x-unknown')],
      ['system', file('image/png')],
      ['assistant', file('image/png')],
    ] as const) {
      await expect(project(adapter, input({ messages: [{ role, parts: [part] }] }))).rejects.toThrow(
        /does not support inline/,
      );
    }
  });

  it('rejects a foreign provider payload with a matching digest', async () => {
    const recorded = recordingFetch(openaiSse('done'));
    const adapter = create(recorded.fetch);
    const request = await project(adapter, input());
    const payload = { ...payloadOf(request), provider: 'google' };
    await expect(collectFrames(adapter, { ...request, payload, digest: await digest(payload) })).rejects.toThrow(
      /provider does not match/,
    );
    expect(recorded.requests).toHaveLength(0);
  });

  it('maps cacheScopeId to promptCacheKey on the wire', async () => {
    const recorded = recordingFetch(openaiSse('done'));
    const adapter = create(recorded.fetch);
    const scoped = input({ configuration: { ...input().configuration, cacheScopeId: 'story-7' } });
    await collectFrames(adapter, await project(adapter, scoped));
    await collectFrames(adapter, await project(adapter, input()));
    expect(JSON.parse(recorded.requests[0]?.body ?? '{}').prompt_cache_key).toBe('story-7');
    expect(JSON.parse(recorded.requests[1]?.body ?? '{}')).not.toHaveProperty('prompt_cache_key');
  });

  it('declares a tool in a system message and calls it in a collected conversation', async () => {
    const recorded = recordingFetch(openaiSse('', { name: 'fixture_lookup', input: '{"id":"x"}' }));
    const adapter = createOpenAIModelAdapter({
      resolveApiKey: () => secret,
      fetch: recorded.fetch,
      modelLimits: { 'gpt-6-custom': { contextWindowTokens: 10_000, maxOutputTokens: 1024 } },
    });
    const request = await project(
      adapter,
      projectInput('openai', 'gpt-6-custom', { messages: [declares(lookupTool), ...input().messages] }),
    );
    expect(resultOf(await collectFrames(adapter, request)).message.parts).toContainEqual(
      expect.objectContaining({ kind: 'tool-call', tool: lookupTool.tool }),
    );
    expect(JSON.parse(recorded.requests[0]?.body ?? '{}').tools[0]).toMatchObject({ name: 'fixture_lookup' });
  });
});
