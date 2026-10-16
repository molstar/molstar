/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type Message, ModelError, type ReasoningParameters } from '../../messages/index.js';
import { createAnthropicProvider } from '../index.js';
import {
  answer,
  anthropicErrorSse,
  anthropicSse,
  anthropicThinkingSse,
  collectFrames,
  hangingFetch,
  lookupTool,
  project,
  projectInput,
  recordingFetch,
  resultOf,
  resultPart,
} from './support.js';

const secret = 'fixture-anthropic-secret';
const create = (fetch?: typeof globalThis.fetch) =>
  createAnthropicProvider({ resolveApiKey: () => secret, ...(fetch === undefined ? {} : { fetch }) });
const sonnet5 = (reasoning?: ReasoningParameters) =>
  projectInput('anthropic', 'claude-sonnet-5', {
    configuration: {
      target: { provider: 'anthropic', model: 'claude-sonnet-5' },
      parameters: { maxOutputTokens: 8000, ...(reasoning === undefined ? {} : { reasoning }) },
    },
  });
const body = (recorded: ReturnType<typeof recordingFetch>, index: number) =>
  JSON.parse(recorded.requests[index]?.body ?? '{}');

describe('Anthropic provider', () => {
  it('streams text and usage, then a call with its signed thinking, which the next request replays first', async () => {
    const text = recordingFetch(anthropicSse('I found it.'));
    const plain = create(text.fetch);
    const frames = await collectFrames(plain, await project(plain, projectInput('anthropic', 'claude-sonnet-4-5')));
    expect(frames).toContainEqual({ kind: 'text-delta', text: 'I found it.' });
    const done = resultOf(frames);
    expect(done).toMatchObject({ finish: 'complete', message: { parts: [{ kind: 'text', text: 'I found it.' }] } });
    // Total input includes the cache reads and writes (12 uncached, 4 read, 1 written).
    expect(done.usage).toMatchObject({ inputTokens: 17, cacheReadInputTokens: 4, cacheWriteInputTokens: 1 });
    expect(text.requests[0]?.headers.get('x-api-key')).toBe(secret);
    expect(text.requests[0]?.body).not.toContain(secret);
    expect(JSON.stringify(frames)).not.toContain(secret);

    const recorded = recordingFetch(
      anthropicThinkingSse(
        { text: 'Look it up.', signature: 'sig-fixture' },
        { name: 'fixture_lookup', input: '{"id":"1abc"}' },
      ),
    );
    const provider = create(recorded.fetch);
    const input = sonnet5({ mode: 'adaptive', effort: 'low' });
    const { message, finish } = resultOf(await collectFrames(provider, await project(provider, input)));
    expect(finish).toBe('tool-calls');
    const call = message.parts.find((part) => part.kind === 'tool-call');
    expect(call).toEqual({
      kind: 'tool-call',
      requestId: 'toolu_fixture',
      tool: lookupTool.tool,
      arguments: { id: '1abc' },
      providerOptions: { anthropic: { thinking: [{ text: 'Look it up.', signature: 'sig-fixture' }] } },
    });
    const next = await project(provider, { ...input, messages: [...input.messages, message, answer('toolu_fixture')] });
    // A retained request (JSON) runs as it is.
    await collectFrames(provider, JSON.parse(JSON.stringify(next)));
    const assistant = body(recorded, 1).messages.find((item: { role: string }) => item.role === 'assistant');
    expect(assistant.content[0]).toEqual({ type: 'thinking', thinking: 'Look it up.', signature: 'sig-fixture' });
    expect(assistant.content.at(-1)).toMatchObject({ type: 'tool_use', id: 'toolu_fixture', name: 'fixture_lookup' });
    expect(body(recorded, 1).messages.at(-1).content[0]).toMatchObject({
      type: 'tool_result',
      tool_use_id: 'toolu_fixture',
    });
  });

  it('replays redacted thinking before its call', async () => {
    const recorded = recordingFetch(
      anthropicThinkingSse({ redactedData: 'opaque-fixture' }, { name: 'fixture_lookup', input: '{"id":"1abc"}' }),
    );
    const provider = create(recorded.fetch);
    const input = sonnet5({ mode: 'adaptive' });
    const { message } = resultOf(await collectFrames(provider, await project(provider, input)));
    expect(message.parts.at(-1)).toMatchObject({
      providerOptions: { anthropic: { thinking: [{ text: '', redactedData: 'opaque-fixture' }] } },
    });
    await collectFrames(
      provider,
      await project(provider, { ...input, messages: [...input.messages, message, answer('toolu_fixture')] }),
    );
    const assistant = body(recorded, 1).messages.find((item: { role: string }) => item.role === 'assistant');
    expect(assistant.content[0]).toEqual({ type: 'redacted_thinking', data: 'opaque-fixture' });
  });

  it('maps a thinking budget and adaptive effort; off is as absent', async () => {
    const recorded = recordingFetch(anthropicSse('Done.'));
    const provider = create(recorded.fetch);
    await collectFrames(provider, await project(provider, sonnet5({ mode: 'budget', tokens: 4_096 })));
    expect(body(recorded, 0).thinking).toMatchObject({ type: 'enabled', budget_tokens: 4_096 });
    await collectFrames(provider, await project(provider, sonnet5({ mode: 'adaptive', effort: 'low' })));
    expect(body(recorded, 1)).toMatchObject({ thinking: { type: 'adaptive' }, output_config: { effort: 'low' } });
    const off = await project(provider, sonnet5({ mode: 'off' }));
    const absent = await project(provider, sonnet5());
    expect(off.payload).toEqual(absent.payload);
  });

  it('refuses reasoning values Claude cannot take, as plain errors', async () => {
    const provider = create();
    for (const reasoning of [
      { mode: 'adaptive', effort: 'minimal' },
      { mode: 'budget', tokens: 512 },
      { mode: 'budget', tokens: 9000 },
    ] as const) {
      const failure = await project(provider, sonnet5(reasoning)).catch((error: unknown) => error);
      expect((failure as Error).message).toMatch(/reasoning/i);
      expect((failure as Error).name).toBe('Error');
    }
  });

  it('puts three cache breakpoints (system, top level, last immediate tool), or none', async () => {
    const guide: Message = { role: 'system', parts: [{ kind: 'text', text: 'Be brief.' }] };
    for (const cache of [true, false]) {
      const recorded = recordingFetch(anthropicSse('ok'));
      const provider = createAnthropicProvider({ resolveApiKey: () => 'fixture', fetch: recorded.fetch, cache });
      const input = projectInput('anthropic', 'claude-sonnet-5-5');
      await collectFrames(provider, await project(provider, { ...input, messages: [guide, ...input.messages] }));
      const wire = body(recorded, 0);
      expect(JSON.stringify(wire).match(/"cache_control"/g)?.length ?? 0).toBe(cache ? 3 : 0);
      if (cache) expect(wire.cache_control).toEqual({ type: 'ephemeral', ttl: '5m' });
    }
  });

  it('sends the tool results first, then the user text around them, in order', async () => {
    const recorded = recordingFetch(anthropicSse('Done.'));
    const provider = create(recorded.fetch);
    const call = (requestId: string) =>
      ({ kind: 'tool-call', requestId, tool: lookupTool.tool, arguments: {} }) as const;
    await collectFrames(
      provider,
      await project(
        provider,
        projectInput('anthropic', 'claude-sonnet-4-5', {
          messages: [
            { role: 'user', parts: [{ kind: 'text', text: 'start' }] },
            { role: 'assistant', parts: [call('toolu_first'), call('toolu_second')] },
            {
              role: 'user',
              parts: [
                { kind: 'text', text: 'before-result' },
                resultPart('toolu_first', 'first result'),
                resultPart('toolu_second', 'second result'),
                { kind: 'text', text: 'after-result' },
              ],
            },
          ],
        }),
      ),
    );
    const wire = recorded.requests[0]?.body ?? '';
    const order = ['first result', 'second result', 'before-result', 'after-result'].map((t) => wire.indexOf(t));
    expect(order.every((at, i) => at > -1 && (i === 0 || at > (order[i - 1] as number)))).toBe(true);
  });

  it('passes its headers, keeps a call to an unknown tool as named, and the finish as the model gave it', async () => {
    const recorded = recordingFetch(
      anthropicSse('', { name: 'Weird_Name-1', input: '{}' }).replace(
        '"stop_reason":"tool_use"',
        '"stop_reason":"end_turn"',
      ),
    );
    const provider = createAnthropicProvider({
      resolveApiKey: () => secret,
      fetch: recorded.fetch,
      headers: { 'x-fixture-header': 'allowed' },
    });
    const result = resultOf(
      await collectFrames(provider, await project(provider, projectInput('anthropic', 'claude-sonnet-4-5'))),
    );
    expect(recorded.requests[0]?.headers.get('x-fixture-header')).toBe('allowed');
    expect(result.finish).toBe('complete');
    expect(result.message.parts).toContainEqual(
      expect.objectContaining({ kind: 'tool-call', tool: { id: 'Weird_Name-1' } }),
    );
  });

  it('aborts an in-flight call, and fails an error event', async () => {
    const hung = hangingFetch();
    const provider = create(hung.fetch);
    const request = await project(provider, projectInput('anthropic', 'claude-sonnet-4-5', { tools: [] }));
    const work = new AbortController();
    const pending = collectFrames(provider, request, work.signal);
    await hung.started;
    work.abort(new Error('test-abort'));
    await expect(pending).rejects.toThrow();
    const failing = create(recordingFetch(anthropicErrorSse()).fetch);
    await expect(collectFrames(failing, request)).rejects.toBeInstanceOf(ModelError);
  });
});
