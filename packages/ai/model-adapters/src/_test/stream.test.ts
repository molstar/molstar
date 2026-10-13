/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ModelError } from '@molstar/ai-agent';
import { MockLanguageModelV4 } from 'ai/test';
import { streamPolicy as anthropicPolicy } from '../anthropic/stream.js';
import { invokeProjectedRequest, mapFinishReason, modelError, providerUsageFromSdk } from '../common/stream.js';
import { streamPolicy as googlePolicy } from '../google/stream.js';
import { createAnthropicModelAdapter, createGoogleModelAdapter, createOpenAIModelAdapter } from '../index.js';
import {
  anthropicSse,
  collectFrames,
  googleSse,
  openaiSse,
  project,
  projectInput,
  recordingFetch,
  resultOf,
} from './support.js';

describe('stream mapping', () => {
  it.each(['anthropic', 'google'] as const)(
    'estimates the input of %s requests from active content',
    async (provider) => {
      const adapter =
        provider === 'anthropic'
          ? createAnthropicModelAdapter({ resolveApiKey: () => 'fixture' })
          : createGoogleModelAdapter({ resolveApiKey: () => 'fixture' });
      const model = provider === 'anthropic' ? 'claude-haiku-4-5' : 'gemini-3.8-flash';
      const text = 'a \\"quoted\\" multilingual é 🙂 context';
      const base = projectInput(provider, model, { messages: [{ role: 'user', parts: [{ kind: 'text', text }] }] });
      const projected = await project(adapter, base);
      expect(projected.inputTokens).not.toBe(JSON.stringify(projected.payload).length);
      const longer = await project(adapter, {
        ...base,
        messages: [...base.messages, { role: 'user', parts: [{ kind: 'text', text: 'append '.repeat(100) }] }],
      });
      expect(longer.inputTokens).toBeGreaterThan(projected.inputTokens + 100);
      // Output and transport controls do not size the prompt.
      const capped = await project(adapter, {
        ...base,
        configuration: { ...base.configuration, parameters: { maxOutputTokens: 32, temperature: 0.5 } },
      });
      expect(capped.inputTokens).toBe(projected.inputTokens);
      expect(capped.digest).not.toBe(projected.digest);
    },
  );

  it.each(['anthropic', 'google'] as const)(
    'preserves an empty %s tool-call finish for consumer validation',
    async (provider) => {
      const adapter =
        provider === 'anthropic'
          ? createAnthropicModelAdapter({ resolveApiKey: () => 'fixture' })
          : createGoogleModelAdapter({ resolveApiKey: () => 'fixture' });
      const request = await project(
        adapter,
        projectInput(provider, provider === 'anthropic' ? 'claude-haiku-4-5' : 'gemini-3.8-flash'),
      );
      const languageModel = new MockLanguageModelV4({
        doStream: {
          stream: new ReadableStream({
            start(controller) {
              controller.enqueue({
                type: 'finish',
                finishReason: { unified: 'tool-calls', raw: 'tool-calls' },
                usage: {
                  inputTokens: { total: 1, noCache: 1, cacheRead: 0, cacheWrite: 0 },
                  outputTokens: { total: 1, text: 1, reasoning: 0 },
                },
              });
              controller.close();
            },
          }),
        },
      });
      const frames = [];
      for await (const frame of invokeProjectedRequest({
        request,
        signal: new AbortController().signal,
        expectedRevision: request.adapterRevision,
        expectedProvider: provider,
        policy: provider === 'anthropic' ? anthropicPolicy : googlePolicy,
        languageModel,
      }))
        frames.push(frame);
      expect(resultOf(frames)).toMatchObject({ finish: 'tool-calls', message: { parts: [] } });
    },
  );

  it('requires a recognized successful finish', () => {
    expect(mapFinishReason('stop')).toBe('complete');
    expect(mapFinishReason('tool-calls')).toBe('tool-calls');
    expect(mapFinishReason('length')).toBe('output-limit');
    expect(mapFinishReason('content-filter')).toBe('filtered');
    expect(() => mapFinishReason('error')).toThrow(/unsuccessfully \(error\)/);
    expect(() => mapFinishReason('other')).toThrow(/unsuccessfully \(other\)/);
    expect(() => mapFinishReason('unknown')).toThrow(/unsuccessfully \(unknown\)/);
    expect(() => mapFinishReason(undefined)).toThrow(/unsuccessfully \(missing\)/);
  });

  it('copies sparse cache and reasoning counts into ProviderUsage', () => {
    expect(
      providerUsageFromSdk({
        inputTokens: 10,
        outputTokens: 2,
        inputTokenDetails: { noCacheTokens: 5, cacheReadTokens: 4, cacheWriteTokens: 1 },
        outputTokenDetails: { textTokens: 1, reasoningTokens: 3 },
        totalTokens: 12,
      }),
    ).toEqual({
      inputTokens: 10,
      outputTokens: 2,
      cacheReadInputTokens: 4,
      cacheWriteInputTokens: 1,
      reasoningOutputTokens: 3,
    });
    expect(
      providerUsageFromSdk({
        inputTokens: undefined,
        outputTokens: undefined,
        inputTokenDetails: { noCacheTokens: undefined, cacheReadTokens: undefined, cacheWriteTokens: undefined },
        outputTokenDetails: { textTokens: undefined, reasoningTokens: undefined },
        totalTokens: undefined,
      }),
    ).toEqual({});
  });

  it.each([
    ['anthropic', createAnthropicModelAdapter, 'claude-sonnet-4-5', anthropicSse],
    ['google', createGoogleModelAdapter, 'gemini-2.5-flash', googleSse],
    ['openai', createOpenAIModelAdapter, 'gpt-6-luna', openaiSse],
  ] as const)(
    'reports %s input as the total with cache reads and writes included',
    async (provider, create, model, sse) => {
      const adapter = create({ resolveApiKey: () => 'fixture', fetch: recordingFetch(sse('done')).fetch });
      const { usage } = resultOf(
        await collectFrames(adapter, await project(adapter, projectInput(provider, model, { tools: [] }))),
      );
      expect(usage.cacheReadInputTokens).toBe(4);
      expect(usage.inputTokens ?? 0).toBeGreaterThanOrEqual(
        (usage.cacheReadInputTokens ?? 0) + (usage.cacheWriteInputTokens ?? 0),
      );
    },
  );
});

describe('failure mapping', () => {
  const fail = (status: number, body: unknown = { error: { message: 'nope' } }) =>
    (async () =>
      new Response(JSON.stringify(body), { status, headers: { 'content-type': 'application/json' } })) as typeof fetch;

  it.each([
    ['anthropic', createAnthropicModelAdapter, 'claude-sonnet-4-5'],
    ['google', createGoogleModelAdapter, 'gemini-2.5-flash'],
    ['openai', createOpenAIModelAdapter, 'gpt-6-luna'],
  ] as const)('maps %s HTTP statuses to transient and effect', async (provider, create, model) => {
    for (const [status, transient] of [
      [400, false],
      [401, false],
      [404, false],
      [408, true],
      [429, true],
      [500, true],
      [503, true],
    ] as const) {
      const adapter = create({ resolveApiKey: () => 'fixture-key', fetch: fail(status) });
      const failure = await collectFrames(adapter, await project(adapter, projectInput(provider, model))).catch(
        (error: unknown) => error,
      );
      expect(failure).toBeInstanceOf(ModelError);
      // The provider refused before generating anything.
      expect(failure).toMatchObject({ transient, effect: 'none' });
    }
  });

  it('treats a network failure as transient with an unknown effect', async () => {
    const adapter = createAnthropicModelAdapter({
      resolveApiKey: () => 'fixture-key',
      fetch: async () => {
        throw new TypeError('fetch failed');
      },
    });
    const failure = await collectFrames(
      adapter,
      await project(adapter, projectInput('anthropic', 'claude-haiku-4-5')),
    ).catch((error: unknown) => error);
    expect(failure).toBeInstanceOf(ModelError);
    expect(failure).toMatchObject({ transient: true, effect: 'unknown' });
  });

  it('refuses a missing API key before any request, without retrying', async () => {
    const adapter = createOpenAIModelAdapter({ resolveApiKey: () => '' });
    const failure = await collectFrames(adapter, await project(adapter, projectInput('openai', 'gpt-6-luna'))).catch(
      (error: unknown) => error,
    );
    expect(failure).toMatchObject({ transient: false, effect: 'none' });
  });

  it('sanitizes the API key out of every failure', async () => {
    const secret = 'sk-very-secret';
    const adapter = createAnthropicModelAdapter({
      resolveApiKey: () => secret,
      fetch: fail(404, { type: 'error', error: { type: 'not_found_error', message: `model ${secret}` } }),
    });
    const failure = (await collectFrames(
      adapter,
      await project(adapter, projectInput('anthropic', 'claude-haiku-4-5')),
    ).catch((error: unknown) => error)) as Error;
    expect(failure.message).toContain('[redacted]');
    expect(failure.message).not.toContain(secret);
    expect(modelError(new Error(`bad ${secret}`), secret).message).toBe('bad [redacted]');
  });

  it("keeps a ModelError's own classification", () => {
    const original = new ModelError('refused', { transient: true, effect: 'none' });
    expect(modelError(original, 'k')).toMatchObject({ message: 'refused', transient: true, effect: 'none' });
    expect(modelError(Object.assign(new Error('x'), { statusCode: 500 }), '')).toMatchObject({
      transient: true,
      effect: 'none',
    });
    expect(modelError(Object.assign(new Error('x'), { isRetryable: false }), '')).toMatchObject({
      transient: false,
      effect: 'unknown',
    });
  });

  it("keeps what the provider said when the SDK did not read its error format (a ChatGPT plan's detail)", () => {
    const body = JSON.stringify({
      detail: "The 'gpt-x' model is not supported when using Codex with a ChatGPT account. k-123",
    });
    const error = modelError(Object.assign(new Error('Bad Request'), { statusCode: 400, responseBody: body }), 'k-123');
    expect(error.message).toBe(
      "Bad Request: The 'gpt-x' model is not supported when using Codex with a ChatGPT account. [redacted] (HTTP 400)",
    );
    expect(error).toMatchObject({ transient: false, effect: 'none' });
    // OpenAI's own format, already in the message, is not repeated.
    const said = JSON.stringify({ error: { message: 'Invalid model' } });
    expect(
      modelError(Object.assign(new Error('Invalid model'), { statusCode: 400, responseBody: said }), '').message,
    ).toBe('Invalid model (HTTP 400)');
  });
});
