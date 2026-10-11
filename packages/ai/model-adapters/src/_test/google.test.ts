/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createGoogleModelAdapter } from '../index.js';
import {
  answer,
  assertDigestMatches,
  collectFrames,
  googleSse,
  hangingFetch,
  lookupTool,
  payloadOf,
  project,
  projectInput,
  recordingFetch,
  resultOf,
  resultPart,
} from './support.js';

const secret = 'fixture-google-secret';
const flash = (extras: Parameters<typeof projectInput>[2] = {}) => projectInput('google', 'gemini-2.5-flash', extras);
const create = (fetch?: typeof globalThis.fetch) =>
  createGoogleModelAdapter({ resolveApiKey: () => secret, ...(fetch === undefined ? {} : { fetch }) });

describe('Google adapter', () => {
  it('maps adaptive effort and token budget to Gemini thinkingConfig', async () => {
    const recorded = recordingFetch(googleSse('Done.'));
    const adapter = create(recorded.fetch);
    const input = projectInput('google', 'gemini-3-flash-preview');
    const withReasoning = (reasoning?: Parameters<typeof project>[1]['configuration']['parameters']['reasoning']) => ({
      ...input,
      configuration: { ...input.configuration, parameters: reasoning === undefined ? {} : { reasoning } },
    });
    const level = await project(adapter, withReasoning({ mode: 'adaptive', effort: 'low' }));
    expect(payloadOf(level).reasoning).toEqual({ mode: 'adaptive', effort: 'low' });
    await assertDigestMatches(level);
    await collectFrames(adapter, level);
    expect(JSON.parse(recorded.requests[0]?.body ?? '{}').generationConfig.thinkingConfig).toEqual({
      thinkingLevel: 'low',
    });
    await collectFrames(adapter, await project(adapter, withReasoning({ mode: 'budget', tokens: 512 })));
    expect(JSON.parse(recorded.requests[1]?.body ?? '{}').generationConfig.thinkingConfig).toEqual({
      thinkingBudget: 512,
    });
    const off = await project(adapter, withReasoning({ mode: 'off' }));
    const absent = await project(adapter, withReasoning());
    expect(off.payload).toEqual(absent.payload);
    expect(off.digest).toBe(absent.digest);
  });

  it('rejects unsupported Gemini reasoning values as plain errors before invocation', async () => {
    const adapter = create();
    const input = projectInput('google', 'gemini-3-flash-preview');
    for (const reasoning of [
      { mode: 'adaptive', effort: 'xhigh' },
      { mode: 'budget', tokens: -1 },
      { mode: 'budget', tokens: 0 },
    ] as const) {
      await expect(
        project(adapter, { ...input, configuration: { ...input.configuration, parameters: { reasoning } } }),
      ).rejects.toThrow(/reasoning/i);
    }
    const older = flash();
    for (const reasoning of [
      { mode: 'adaptive', effort: 'low' },
      { mode: 'budget', tokens: 512 },
    ] as const) {
      await expect(
        project(adapter, { ...older, configuration: { ...older.configuration, parameters: { reasoning } } }),
      ).rejects.toThrow(/Gemini 3 or later/);
    }
  });

  it('advertises the output limit used when none is configured, and no tool appending or cache scope', async () => {
    const adapter = create();
    const configuration = { ...flash().configuration, parameters: {} };
    const capabilities = await adapter.capabilitiesFor(configuration);
    expect(payloadOf(await project(adapter, { ...flash(), configuration })).maxOutputTokens).toBe(4096);
    expect(capabilities).toMatchObject({
      defaultMaxOutputTokens: 4096,
      maxInputTokens: 1_048_576,
      appendsTools: false,
      cacheScope: false,
      outputLimitIncludesReasoning: false,
    });
  });

  it('sends all mixed user text after the required adjacent tool result', async () => {
    const recorded = recordingFetch(googleSse('Done.'));
    const adapter = create(recorded.fetch);
    const projected = await project(
      adapter,
      flash({
        messages: [
          { role: 'user', parts: [{ kind: 'text', text: 'start' }] },
          {
            role: 'assistant',
            parts: [{ kind: 'tool-call', requestId: 'call_fixture', tool: lookupTool.tool, arguments: {} }],
          },
          {
            role: 'user',
            parts: [
              { kind: 'text', text: 'before-result' },
              resultPart('call_fixture'),
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

  it('projects a stable digest without placing the API key in the payload', async () => {
    const adapter = create();
    const first = await project(adapter, flash());
    const second = await project(adapter, flash());
    await assertDigestMatches(first);
    expect(first.digest).toBe(second.digest);
    expect(first.adapterRevision).toBe('google-ai-sdk-v3');
    expect(JSON.stringify(first.payload)).not.toContain(secret);
  });

  it('streams text and a final result without leaking the API key into retained frames', async () => {
    const recorded = recordingFetch(googleSse('Ordinary response.'));
    const adapter = create(recorded.fetch);
    const frames = await collectFrames(adapter, await project(adapter, flash({ tools: [] })));
    expect(frames.some((frame) => frame.kind === 'text-delta' && frame.text.includes('Ordinary response.'))).toBe(true);
    const result = resultOf(frames);
    expect(result.finish).toBe('complete');
    // Gemini's prompt count already includes the cached tokens.
    expect(result.usage).toMatchObject({ inputTokens: 12, cacheReadInputTokens: 4 });
    expect(JSON.stringify(frames)).not.toContain(secret);
    expect(recorded.requests[0]?.headers.get('x-goog-api-key')).toBe(secret);
    expect(recorded.requests[0]?.body).not.toContain(secret);
  });

  it('maps streamed function calls back to the original tool id', async () => {
    const recorded = recordingFetch(googleSse('I found it.', { name: 'fixture_lookup', args: { id: '1abc' } }));
    const adapter = create(recorded.fetch);
    const result = resultOf(await collectFrames(adapter, await project(adapter, flash())));
    expect(result.finish).toBe('tool-calls');
    expect(result.message.parts).toContainEqual({
      kind: 'tool-call',
      requestId: 'call_fixture',
      tool: { id: lookupTool.tool.id },
      arguments: { id: '1abc' },
    });
  });

  it('retains a function-call thought signature through the next Google request', async () => {
    const firstFetch = recordingFetch(
      googleSse('Looking up the structure.', {
        name: 'fixture_lookup',
        args: { id: '1abc' },
        thoughtSignature: 'signed_fixture',
      }),
    );
    const first = create(firstFetch.fetch);
    const model = 'gemini-3-flash-preview';
    const { message } = resultOf(await collectFrames(first, await project(first, projectInput('google', model))));
    expect(message.parts).toContainEqual({
      kind: 'tool-call',
      requestId: 'call_fixture',
      tool: { id: lookupTool.tool.id },
      arguments: { id: '1abc' },
      providerOptions: { google: { thoughtSignature: 'signed_fixture' } },
    });
    const followUpFetch = recordingFetch(googleSse('The structure is available.'));
    const followUpAdapter = create(followUpFetch.fetch);
    const followUp = await project(
      followUpAdapter,
      projectInput('google', model, {
        messages: [
          { role: 'user', parts: [{ kind: 'text', text: 'Look up 1abc.' }] },
          message,
          answer('call_fixture', 'Found 1abc.'),
        ],
      }),
    );
    await assertDigestMatches(followUp);
    expect(JSON.stringify(followUp.payload)).toContain('signed_fixture');
    await collectFrames(followUpAdapter, followUp);
    expect(followUpFetch.requests[0]?.body).toContain('"thoughtSignature":"signed_fixture"');
    expect(followUpFetch.requests[0]?.body).not.toContain('skip_thought_signature_validator');
  });

  it('aborts an in-flight invoke', async () => {
    const hung = hangingFetch();
    const adapter = create(hung.fetch);
    const projected = await project(adapter, flash({ tools: [] }));
    const work = new AbortController();
    const pending = collectFrames(adapter, projected, work.signal);
    await hung.started;
    work.abort(new Error('test-abort'));
    await expect(pending).rejects.toThrow();
  });
});
