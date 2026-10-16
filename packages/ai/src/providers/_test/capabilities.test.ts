/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createAnthropicProvider, createGoogleProvider, createOpenAIProvider } from '../index.js';
import { anthropicSse, collectFrames, payloadOf, project, projectInput, recordingFetch } from './support.js';

describe('effective model capabilities', () => {
  it('resolves limits for different models of one provider', async () => {
    const provider = createAnthropicProvider({ resolveApiKey: () => 'unused' });
    const small = await provider.capabilities(projectInput('anthropic', 'claude-haiku-4-5').configuration);
    const large = await provider.capabilities(projectInput('anthropic', 'claude-sonnet-5').configuration);
    expect(small.contextWindowTokens).toBe(200_000);
    expect(large.contextWindowTokens).toBe(1_000_000);
    expect(small.maxOutputTokens).toBe(64_000);
    expect(large.maxOutputTokens).toBe(128_000);
  });

  it.each([
    ['anthropic', 'claude-sonnet-5-5', createAnthropicProvider, { appendsTools: true, cacheScope: false }],
    ['anthropic', 'claude-sonnet-5', createAnthropicProvider, { appendsTools: false, cacheScope: false }],
    ['openai', 'gpt-6-luna', createOpenAIProvider, { appendsTools: true, cacheScope: true }],
    ['google', 'gemini-3.8-flash', createGoogleProvider, { appendsTools: false, cacheScope: false }],
  ] as const)('reports what %s %s can do', async (vendor, model, create, expected) => {
    const provider = create({ resolveApiKey: () => 'unused' });
    expect(await provider.capabilities(projectInput(vendor, model).configuration)).toMatchObject(expected);
  });

  it('rejects unknown targets and allows explicit gateway limits', async () => {
    const configuration = { target: { provider: 'google', model: 'custom-gateway' }, parameters: {} };
    const unknown = createGoogleProvider({ resolveApiKey: () => 'unused' });
    await expect(unknown.capabilities(configuration)).rejects.toThrow('explicit token limits');
    const custom = createGoogleProvider({
      resolveApiKey: () => 'unused',
      modelLimits: {
        'custom-gateway': { contextWindowTokens: 20_000, maxOutputTokens: 2048, defaultMaxOutputTokens: 512 },
      },
    });
    const capabilities = await custom.capabilities(configuration);
    expect(capabilities.contextWindowTokens).toBe(20_000);
    const request = await custom.project({ configuration, messages: [], tools: [], capabilities });
    expect(payloadOf(request).maxOutputTokens).toBe(512);
    await expect(
      custom.project({
        configuration: { ...configuration, parameters: { maxOutputTokens: 2049 } },
        messages: [],
        tools: [],
        capabilities,
      }),
    ).rejects.toThrow('exceeds model maximum');
  });

  it("projects against the capabilities it is given, not the factory's current overrides", async () => {
    const input = projectInput('google', 'gemini-3.8-flash');
    const first = createGoogleProvider({
      resolveApiKey: () => 'unused',
      modelLimits: { 'gemini-3.8-flash': { contextWindowTokens: 100_000, maxOutputTokens: 4000 } },
    });
    const capabilities = await first.capabilities(input.configuration);
    const original = await first.project({ ...input, capabilities });
    const second = createGoogleProvider({
      resolveApiKey: () => 'unused',
      modelLimits: { 'gemini-3.8-flash': { contextWindowTokens: 10_000, maxOutputTokens: 32 } },
    });
    expect(await second.project({ ...input, capabilities })).toEqual(original);
  });

  it('keeps prompt sizing independent of output and transport controls', async () => {
    const provider = createGoogleProvider({ resolveApiKey: () => 'unused' });
    const input = projectInput('google', 'gemini-3.8-flash');
    const first = await project(provider, input);
    const second = await project(provider, {
      ...input,
      configuration: { ...input.configuration, parameters: { maxOutputTokens: 32, temperature: 0.5 } },
    });
    expect(first.inputTokens).toBe(second.inputTokens);
    expect(first.digest).not.toBe(second.digest);
  });

  it('rejects a default output limit above the model maximum', async () => {
    const provider = createGoogleProvider({
      resolveApiKey: () => 'unused',
      modelLimits: { odd: { contextWindowTokens: 1000, maxOutputTokens: 100, defaultMaxOutputTokens: 200 } },
    });
    await expect(
      provider.capabilities({ target: { provider: 'google', model: 'odd' }, parameters: {} }),
    ).rejects.toThrow('Default output limit');
  });
});

// Published count_tokens outputs, captured 2026-09-30 from the primary examples:
// https://platform.claude.com/docs/en/build-with-claude/token-counting
// These are calibration samples, not synthetic SSE usage or exact-tokenizer claims.
describe('published token-count calibration', () => {
  it('stays within a broad error band for prose and a tool schema', async () => {
    const provider = createAnthropicProvider({ resolveApiKey: () => 'unused' });
    const configuration = { target: { provider: 'anthropic', model: 'claude-sonnet-4-5' }, parameters: {} };
    const prose = await project(provider, {
      configuration,
      messages: [
        { role: 'system', parts: [{ kind: 'text', text: 'You are a scientist' }] },
        { role: 'user', parts: [{ kind: 'text', text: 'Hello, Claude' }] },
      ],
      tools: [],
    });
    const schema = await project(provider, {
      configuration,
      messages: [{ role: 'user', parts: [{ kind: 'text', text: "What's the weather like in San Francisco?" }] }],
      tools: [
        {
          tool: { id: 'get_weather' },
          name: 'get_weather',
          description: 'Get the current weather in a given location',
          inputSchema: {
            type: 'object',
            properties: { location: { type: 'string', description: 'The city and state, e.g. San Francisco, CA' } },
            required: ['location'],
          },
        },
      ],
    });
    for (const [estimate, recorded] of [
      [prose.inputTokens, 14],
      [schema.inputTokens, 403],
    ] as const) {
      expect(estimate / recorded).toBeGreaterThan(0.7);
      expect(estimate / recorded).toBeLessThan(1.5);
    }
  });

  it('is not needed by an invocation: a retained request runs on its own payload', async () => {
    const recorded = recordingFetch(anthropicSse('ok'));
    const provider = createAnthropicProvider({ resolveApiKey: () => 'k', fetch: recorded.fetch });
    const request = await project(provider, projectInput('anthropic', 'claude-haiku-4-5'));
    await collectFrames(provider, JSON.parse(JSON.stringify(request)));
    expect(recorded.requests).toHaveLength(1);
  });
});
