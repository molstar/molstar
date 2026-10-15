/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ModelError, type ReasoningParameters } from '@molstar/ai-agent';
import { createAnthropicModelAdapter } from '../index.js';
import {
  answer,
  anthropicErrorSse,
  anthropicSse,
  anthropicThinkingSse,
  assertDigestMatches,
  collectFrames,
  hangingFetch,
  lookupTool,
  payloadOf,
  project,
  projectInput,
  recordingFetch,
  resultOf,
  resultPart,
} from './support.js';

const secret = 'fixture-anthropic-secret';
const sonnet45 = (extras: Parameters<typeof projectInput>[2] = {}) =>
  projectInput('anthropic', 'claude-sonnet-4-5', extras);
const create = (fetch?: typeof globalThis.fetch) =>
  createAnthropicModelAdapter({ resolveApiKey: () => secret, ...(fetch === undefined ? {} : { fetch }) });

describe('Anthropic adapter', () => {
  it('projects the selected reasoning budget into the provider request', async () => {
    const recorded = recordingFetch(anthropicSse('Done.'));
    const adapter = create(recorded.fetch);
    const input = sonnet45();
    const projected = await project(adapter, {
      ...input,
      configuration: {
        ...input.configuration,
        parameters: { maxOutputTokens: 8_000, reasoning: { mode: 'budget', tokens: 4_096 } },
      },
    });
    expect(payloadOf(projected).reasoning).toEqual({ mode: 'budget', tokens: 4_096 });
    await assertDigestMatches(projected);
    await collectFrames(adapter, projected);
    const body = JSON.parse(recorded.requests[0]?.body ?? '{}');
    expect(body.thinking).toMatchObject({ type: 'enabled', budget_tokens: 4_096 });
  });

  it('advertises the output limit used when none is configured', async () => {
    const adapter = create();
    const input = sonnet45();
    const configuration = { ...input.configuration, parameters: {} };
    const capabilities = await adapter.capabilitiesFor(configuration);
    const projected = await project(adapter, { ...input, configuration });
    expect(payloadOf(projected).maxOutputTokens).toBe(4096);
    expect(capabilities.defaultMaxOutputTokens).toBe(4096);
    expect(capabilities.maxOutputTokens).toBeGreaterThan(4096);
  });

  it('sends all mixed user text after the required adjacent tool result', async () => {
    const recorded = recordingFetch(anthropicSse('Done.'));
    const adapter = create(recorded.fetch);
    const call = { kind: 'tool-call', requestId: 'toolu_fixture', tool: lookupTool.tool, arguments: {} } as const;
    const projected = await project(
      adapter,
      sonnet45({
        messages: [
          { role: 'user', parts: [{ kind: 'text', text: 'start' }] },
          { role: 'assistant', parts: [call] },
          {
            role: 'user',
            parts: [
              { kind: 'text', text: 'before-result' },
              resultPart('toolu_fixture'),
              { kind: 'text', text: 'after-result' },
            ],
          },
        ],
      }),
    );
    expect((payloadOf(projected).messages as readonly { role: string }[]).map((message) => message.role)).toEqual([
      'user',
      'assistant',
      'tool',
      'user',
    ]);
    await collectFrames(adapter, projected);
    const body = recorded.requests[0]?.body ?? '';
    expect(body.indexOf('found')).toBeLessThan(body.indexOf('before-result'));
    expect(body.indexOf('before-result')).toBeLessThan(body.indexOf('after-result'));
  });

  it('sends adjacent tool results in the order supplied', async () => {
    const recorded = recordingFetch(anthropicSse('Done.'));
    const adapter = create(recorded.fetch);
    const projected = await project(
      adapter,
      sonnet45({
        messages: [
          {
            role: 'assistant',
            parts: [
              { kind: 'tool-call', requestId: 'toolu_first', tool: lookupTool.tool, arguments: {} },
              { kind: 'tool-call', requestId: 'toolu_second', tool: lookupTool.tool, arguments: {} },
            ],
          },
          {
            role: 'user',
            parts: [resultPart('toolu_first', 'first result'), resultPart('toolu_second', 'second result')],
          },
        ],
      }),
    );
    await collectFrames(adapter, projected);
    const body = recorded.requests[0]?.body ?? '';
    expect(body.indexOf('first result')).toBeGreaterThan(-1);
    expect(body.indexOf('first result')).toBeLessThan(body.indexOf('second result'));
  });

  it('projects a stable digest without placing the API key in the payload', async () => {
    const adapter = create();
    const first = await project(adapter, sonnet45());
    const second = await project(adapter, sonnet45());
    await assertDigestMatches(first);
    expect(first.digest).toBe(second.digest);
    expect(first.digest).toMatch(/^[a-f0-9]{64}$/);
    expect(first.adapterRevision).toBe('anthropic-ai-sdk-v3');
    expect(JSON.stringify(first.payload)).not.toContain(secret);
    expect(payloadOf(first).model).toBe('claude-sonnet-4-5');
    expect(JSON.stringify(payloadOf(first).tools)).toContain('fixture_lookup');
    expect(first.inputTokens).toBeGreaterThan(0);
  });

  it('maps adaptive thinking and effort into the Claude request', async () => {
    const recorded = recordingFetch(anthropicSse('Done.'));
    const adapter = create(recorded.fetch);
    const input = projectInput('anthropic', 'claude-sonnet-5');
    const withReasoning = (reasoning: ReasoningParameters) => ({
      ...input,
      configuration: { ...input.configuration, parameters: { maxOutputTokens: 8000, reasoning } },
    });
    const adaptive = await project(adapter, withReasoning({ mode: 'adaptive', effort: 'low' }));
    expect(payloadOf(adaptive).reasoning).toEqual({ mode: 'adaptive', effort: 'low' });
    await collectFrames(adapter, adaptive);
    expect(JSON.parse(recorded.requests[0]?.body ?? '{}')).toMatchObject({
      thinking: { type: 'adaptive' },
      output_config: { effort: 'low' },
    });
    const off = await project(adapter, withReasoning({ mode: 'off' }));
    const absent = await project(adapter, {
      ...input,
      configuration: { ...input.configuration, parameters: { maxOutputTokens: 8000 } },
    });
    expect(off.payload).toEqual(absent.payload);
    expect(off.digest).toBe(absent.digest);
  });

  it('rejects unsupported Claude reasoning values and a temperature-free parameter set as plain errors', async () => {
    const adapter = create();
    const input = projectInput('anthropic', 'claude-sonnet-5');
    for (const reasoning of [
      { mode: 'adaptive', effort: 'minimal' },
      { mode: 'budget', tokens: 512 },
    ] as const) {
      const failure = await project(adapter, {
        ...input,
        configuration: { ...input.configuration, parameters: { reasoning } },
      }).catch((error: unknown) => error);
      expect(failure).toBeInstanceOf(Error);
      expect((failure as Error).message).toMatch(/reasoning/i);
      expect((failure as Error).name).toBe('Error');
    }
  });

  it('keeps signed Claude thinking on the following tool call and replays it before the call', async () => {
    const recorded = recordingFetch(
      anthropicThinkingSse(
        { text: 'Look it up.', signature: 'sig-fixture' },
        { name: 'fixture_lookup', input: '{"id":"1abc"}' },
      ),
    );
    const adapter = create(recorded.fetch);
    const input = projectInput('anthropic', 'claude-sonnet-5', {
      configuration: {
        target: { provider: 'anthropic', model: 'claude-sonnet-5' },
        parameters: { maxOutputTokens: 8000, reasoning: { mode: 'adaptive', effort: 'low' } },
      },
    });
    const { message } = resultOf(await collectFrames(adapter, await project(adapter, input)));
    const call = message.parts.find((part) => part.kind === 'tool-call');
    if (call?.kind !== 'tool-call') throw new Error('Expected a tool call');
    expect(call.providerOptions).toEqual({
      anthropic: { thinking: [{ text: 'Look it up.', signature: 'sig-fixture' }] },
    });
    const next = await project(adapter, { ...input, messages: [...input.messages, message, answer(call.requestId)] });
    await assertDigestMatches(next);
    await collectFrames(adapter, next);
    const body = JSON.parse(recorded.requests[1]?.body ?? '{}');
    const assistant = body.messages.find((item: { role: string }) => item.role === 'assistant');
    expect(assistant.content[0]).toEqual({ type: 'thinking', thinking: 'Look it up.', signature: 'sig-fixture' });
    expect(assistant.content.at(-1).type).toBe('tool_use');
  });

  it('replays redacted Claude thinking before its tool call', async () => {
    const recorded = recordingFetch(
      anthropicThinkingSse({ redactedData: 'opaque-fixture' }, { name: 'fixture_lookup', input: '{"id":"1abc"}' }),
    );
    const adapter = create(recorded.fetch);
    const input = projectInput('anthropic', 'claude-sonnet-5', {
      configuration: {
        target: { provider: 'anthropic', model: 'claude-sonnet-5' },
        parameters: { maxOutputTokens: 8000, reasoning: { mode: 'adaptive' } },
      },
    });
    const { message } = resultOf(await collectFrames(adapter, await project(adapter, input)));
    const call = message.parts.find((part) => part.kind === 'tool-call');
    if (call?.kind !== 'tool-call') throw new Error('Expected a tool call');
    expect(call.providerOptions).toEqual({ anthropic: { thinking: [{ text: '', redactedData: 'opaque-fixture' }] } });
    await collectFrames(
      adapter,
      await project(adapter, { ...input, messages: [...input.messages, message, answer(call.requestId)] }),
    );
    const body = JSON.parse(recorded.requests[1]?.body ?? '{}');
    const assistant = body.messages.find((item: { role: string }) => item.role === 'assistant');
    expect(assistant.content[0]).toEqual({ type: 'redacted_thinking', data: 'opaque-fixture' });
    expect(assistant.content.at(-1).type).toBe('tool_use');
  });

  it('passes optional Anthropic headers to the provider', async () => {
    const recorded = recordingFetch(anthropicSse('Done.'));
    const adapter = createAnthropicModelAdapter({
      resolveApiKey: () => secret,
      fetch: recorded.fetch,
      headers: { 'x-fixture-header': 'allowed' },
    });
    await collectFrames(adapter, await project(adapter, sonnet45()));
    expect(recorded.requests[0]?.headers.get('x-fixture-header')).toBe('allowed');
  });

  it('derives the callable alias from id rather than display name', async () => {
    const adapter = create();
    const projected = await project(adapter, sonnet45({ tools: [{ ...lookupTool, name: 'Pretty label' }] }));
    expect(payloadOf(projected).tools).toMatchObject([{ id: 'fixture.lookup', wireName: 'fixture_lookup' }]);
    expect(JSON.stringify(projected.payload)).not.toContain('Pretty label');
  });

  it('streams text, usage with total input, and a final result without leaking the API key', async () => {
    const recorded = recordingFetch(anthropicSse('I found it.'));
    const adapter = create(recorded.fetch);
    const frames = await collectFrames(adapter, await project(adapter, sonnet45({ tools: [] })));
    expect(frames.some((frame) => frame.kind === 'text-delta' && frame.text.includes('I found it.'))).toBe(true);
    const result = resultOf(frames);
    expect(result.finish).toBe('complete');
    expect(result.message.parts).toContainEqual({ kind: 'text', text: 'I found it.' });
    expect(result.usage.outputTokens).toBeGreaterThan(0);
    // Total input includes the cache reads and writes (12 uncached, 4 read, 1 written).
    expect(result.usage.inputTokens).toBe(17);
    expect(result.usage).toMatchObject({ cacheReadInputTokens: 4, cacheWriteInputTokens: 1 });
    expect(JSON.stringify(frames)).not.toContain(secret);
    expect(recorded.requests).toHaveLength(1);
    expect(recorded.requests[0]?.headers.get('x-api-key')).toBe(secret);
    expect(recorded.requests[0]?.body).not.toContain(secret);
  });

  it('maps streamed tool calls back to the original tool id', async () => {
    const recorded = recordingFetch(anthropicSse('I found it.', { name: 'fixture_lookup', input: '{"id":"1abc"}' }));
    const adapter = create(recorded.fetch);
    const result = resultOf(await collectFrames(adapter, await project(adapter, sonnet45())));
    expect(result.finish).toBe('tool-calls');
    expect(result.message.parts).toContainEqual({
      kind: 'tool-call',
      requestId: 'toolu_fixture',
      tool: { id: lookupTool.tool.id },
      arguments: { id: '1abc' },
    });
    expect(JSON.parse(recorded.requests[0]?.body ?? '{}').tools[0].name).toBe('fixture_lookup');
  });

  it('keeps the name of a call to an unknown tool as the model gave it', async () => {
    const recorded = recordingFetch(anthropicSse('', { name: 'Weird_Name-1', input: '{}' }));
    const adapter = create(recorded.fetch);
    const result = resultOf(await collectFrames(adapter, await project(adapter, sonnet45())));
    expect(result.message.parts).toContainEqual(
      expect.objectContaining({ kind: 'tool-call', tool: { id: 'Weird_Name-1' } }),
    );
  });

  it('preserves a tool-use finish without a collected tool call for consumer validation', async () => {
    const recorded = recordingFetch(
      anthropicSse('No call.').replace('"stop_reason":"end_turn"', '"stop_reason":"tool_use"'),
    );
    const adapter = create(recorded.fetch);
    const result = resultOf(await collectFrames(adapter, await project(adapter, sonnet45())));
    expect(result).toMatchObject({ finish: 'tool-calls', message: { parts: [{ kind: 'text' }] } });
  });

  it('preserves a complete finish with a collected tool call for consumer validation', async () => {
    const recorded = recordingFetch(
      anthropicSse('I found it.', { name: 'fixture_lookup', input: '{"id":"1abc"}' }).replace(
        '"stop_reason":"tool_use"',
        '"stop_reason":"end_turn"',
      ),
    );
    const adapter = create(recorded.fetch);
    const result = resultOf(await collectFrames(adapter, await project(adapter, sonnet45())));
    expect(result.finish).toBe('complete');
    expect(result.message.parts.some((part) => part.kind === 'tool-call')).toBe(true);
  });

  it('sends tool results with the original call id on the next invoke', async () => {
    const recorded = recordingFetch(anthropicSse('Done.'));
    const adapter = create(recorded.fetch);
    const projected = await project(
      adapter,
      sonnet45({
        messages: [
          { role: 'user', parts: [{ kind: 'text', text: 'Look up 1abc.' }] },
          {
            role: 'assistant',
            parts: [
              { kind: 'text', text: 'calling' },
              { kind: 'tool-call', requestId: 'toolu_fixture', tool: lookupTool.tool, arguments: { id: '1abc' } },
            ],
          },
          answer('toolu_fixture'),
        ],
      }),
    );
    expect(resultOf(await collectFrames(adapter, projected)).finish).toBe('complete');
    expect(JSON.stringify(projected.payload)).toContain('toolu_fixture');
    expect(JSON.stringify(payloadOf(projected).messages)).toContain('"role":"tool"');
  });

  it('rejects a tampered payload digest, a foreign revision and a foreign provider before any request', async () => {
    const recorded = recordingFetch(anthropicSse('Done.'));
    const adapter = create(recorded.fetch);
    const projected = await project(adapter, sonnet45({ tools: [] }));
    const refusals = [
      [{ ...projected, digest: 'tampered' }, /payload digest/i],
      [{ ...projected, adapterRevision: 'old' }, /incompatible adapter revision/i],
    ] as const;
    for (const [request, message] of refusals) {
      const failure = await collectFrames(adapter, request).catch((error: unknown) => error);
      expect(failure).toBeInstanceOf(ModelError);
      expect(failure).toMatchObject({ transient: false, effect: 'none' });
      expect((failure as Error).message).toMatch(message);
    }
    expect(recorded.requests).toHaveLength(0);
  });

  it('aborts an in-flight invoke', async () => {
    const hung = hangingFetch();
    const adapter = create(hung.fetch);
    const projected = await project(adapter, sonnet45({ tools: [] }));
    const work = new AbortController();
    const pending = collectFrames(adapter, projected, work.signal);
    await hung.started;
    work.abort(new Error('test-abort'));
    await expect(pending).rejects.toThrow();
  });

  it('fails an error event instead of accepting a complete result', async () => {
    const adapter = create(recordingFetch(anthropicErrorSse()).fetch);
    await expect(collectFrames(adapter, await project(adapter, sonnet45({ tools: [] })))).rejects.toBeInstanceOf(
      ModelError,
    );
  });

  it('rejects a configuration for another provider and an unknown model', async () => {
    const adapter = create();
    await expect(project(adapter, projectInput('openai' as never, 'gpt-6-luna'))).rejects.toThrow(/cannot project/);
    await expect(
      adapter.capabilitiesFor(
        sonnet45({ configuration: { target: { provider: 'anthropic', model: 'x' }, parameters: {} } }).configuration,
      ),
    ).rejects.toThrow(/explicit token limits/);
  });
});
