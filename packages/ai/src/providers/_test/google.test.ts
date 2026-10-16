/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { ReasoningParameters } from '../../messages/index.js';
import { createGoogleProvider } from '../index.js';
import {
  answer,
  collectFrames,
  googleSse,
  lookupTool,
  payloadOf,
  project,
  projectInput,
  recordingFetch,
  resultOf,
} from './support.js';

const secret = 'fixture-google-secret';
const create = (fetch?: typeof globalThis.fetch) =>
  createGoogleProvider({ resolveApiKey: () => secret, ...(fetch === undefined ? {} : { fetch }) });
const gemini3 = (reasoning?: ReasoningParameters, model = 'gemini-3-flash-preview') =>
  projectInput('google', model, {
    configuration: {
      target: { provider: 'google', model },
      parameters: reasoning === undefined ? {} : { reasoning },
    },
  });

describe('Google provider', () => {
  it('streams text and usage, then a call with its thought signature, which the next request sends back', async () => {
    const text = recordingFetch(googleSse('Ordinary response.'));
    const plain = create(text.fetch);
    const frames = await collectFrames(plain, await project(plain, projectInput('google', 'gemini-2.5-flash')));
    expect(frames).toContainEqual({ kind: 'text-delta', text: 'Ordinary response.' });
    // Gemini's prompt count already includes the cached tokens.
    expect(resultOf(frames)).toMatchObject({ finish: 'complete', usage: { inputTokens: 12, cacheReadInputTokens: 4 } });
    expect(text.requests[0]?.headers.get('x-goog-api-key')).toBe(secret);
    expect(text.requests[0]?.body).not.toContain(secret);

    const recorded = recordingFetch(
      googleSse('Looking it up.', { name: 'fixture_lookup', args: { id: '1abc' }, thoughtSignature: 'signed_fixture' }),
    );
    const provider = create(recorded.fetch);
    const input = gemini3();
    const { message, finish } = resultOf(await collectFrames(provider, await project(provider, input)));
    expect(finish).toBe('tool-calls');
    expect(message.parts).toContainEqual({
      kind: 'tool-call',
      requestId: 'call_fixture',
      tool: lookupTool.tool,
      arguments: { id: '1abc' },
      providerOptions: { google: { thoughtSignature: 'signed_fixture' } },
    });
    await collectFrames(
      provider,
      await project(provider, { ...input, messages: [...input.messages, message, answer('call_fixture', 'Found.')] }),
    );
    const sent = JSON.parse(recorded.requests[1]?.body ?? '{}');
    expect(sent.contents[1].parts.at(-1)).toMatchObject({
      functionCall: { name: 'fixture_lookup', args: { id: '1abc' } },
      thoughtSignature: 'signed_fixture',
    });
    expect(sent.contents[2].parts[0].functionResponse).toMatchObject({ name: 'fixture_lookup' });
    expect(recorded.requests[1]?.body).not.toContain('skip_thought_signature_validator');
  });

  it('maps an effort and a budget to thinkingConfig; off is as absent', async () => {
    const recorded = recordingFetch(googleSse('Done.'));
    const provider = create(recorded.fetch);
    await collectFrames(provider, await project(provider, gemini3({ mode: 'adaptive', effort: 'low' })));
    await collectFrames(provider, await project(provider, gemini3({ mode: 'budget', tokens: 512 })));
    const config = (i: number) => JSON.parse(recorded.requests[i]?.body ?? '{}').generationConfig.thinkingConfig;
    expect(config(0)).toEqual({ thinkingLevel: 'low' });
    expect(config(1)).toEqual({ thinkingBudget: 512 });
    expect((await project(provider, gemini3({ mode: 'off' }))).payload).toEqual(
      (await project(provider, gemini3())).payload,
    );
  });

  it('refuses reasoning values Gemini cannot take, and any before Gemini 3, as plain errors', async () => {
    const provider = create();
    for (const reasoning of [
      { mode: 'adaptive', effort: 'xhigh' },
      { mode: 'budget', tokens: 0 },
    ] as const) {
      await expect(project(provider, gemini3(reasoning))).rejects.toThrow(/reasoning/i);
    }
    await expect(project(provider, gemini3({ mode: 'adaptive', effort: 'low' }, 'gemini-2.5-flash'))).rejects.toThrow(
      /Gemini 3 or later/,
    );
  });

  it('collects its tools, has no cache scope, and sends the default output limit', async () => {
    const provider = create();
    const input = projectInput('google', 'gemini-2.5-flash');
    const configuration = { ...input.configuration, parameters: {} };
    expect(await provider.capabilities(configuration)).toMatchObject({
      defaultMaxOutputTokens: 4096,
      maxInputTokens: 1_048_576,
      appendsTools: false,
      cacheScope: false,
      outputLimitIncludesReasoning: false,
    });
    expect(payloadOf(await project(provider, { ...input, configuration })).maxOutputTokens).toBe(4096);
  });
});
