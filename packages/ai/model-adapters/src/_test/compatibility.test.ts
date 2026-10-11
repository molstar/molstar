/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import type { Message, ModelCapabilities, ProjectInput } from '@molstar/ai-agent';
import { createAnthropicModelAdapter, createGoogleModelAdapter, createOpenAIModelAdapter } from '../index.js';
import {
  anthropicSse,
  collectFrames,
  declares,
  googleSse,
  lookupTool,
  openaiSse,
  type Provider,
  recordingFetch,
  resultOf,
} from './support.js';

interface WireFixture {
  readonly name: string;
  readonly provider: Provider;
  readonly mode: 'collected' | 'reasoning' | 'append-only';
  readonly wireDigests: readonly string[];
}
// HTTP body fingerprints captured from the vendored adapters before the port to the new agent
// package: the same conversation must still produce the same bytes on the wire.
const fixtures = JSON.parse(readFileSync(join(__dirname, 'retained-wire.json'), 'utf8')) as {
  readonly fixtures: readonly WireFixture[];
};
const factories = {
  anthropic: createAnthropicModelAdapter,
  google: createGoogleModelAdapter,
  openai: createOpenAIModelAdapter,
};
const responses = { anthropic: anthropicSse, google: googleSse, openai: openaiSse };
const models = { anthropic: 'claude-sonnet-5-5', google: 'gemini-3.8-flash', openai: 'gpt-6-luna' } as const;
const second = { ...lookupTool, tool: { id: 'second-tool' }, name: 'second-tool' };

function conversation(): Message[] {
  return [
    {
      role: 'system',
      parts: [
        { kind: 'tool-definition', definition: lookupTool, deferred: false },
        { kind: 'text', text: 'Initial guidance' },
      ],
    },
    {
      role: 'user',
      parts: [
        { kind: 'text', text: 'Read é 🙂' },
        { kind: 'reference', reference: { uri: 'data:image/png;base64,Ynl0ZXM=' } },
        { kind: 'reference', reference: { uri: 'data:text/plain;base64,dGV4dA==' } },
      ],
    },
    {
      role: 'assistant',
      parts: [
        {
          kind: 'tool-call',
          tool: lookupTool.tool,
          requestId: 'call-fixture',
          arguments: { id: 'fixture' },
          providerOptions: {
            anthropic: { thinking: [{ text: 'reason', signature: 'signature' }] },
            google: { thoughtSignature: 'signature' },
            openai: {
              itemId: 'fc-fixture',
              reasoning: [{ itemId: 'rs-fixture', text: 'reason', encryptedContent: 'encrypted' }],
            },
          },
        },
      ],
    },
    {
      role: 'user',
      parts: [
        {
          kind: 'tool-result',
          tool: lookupTool.tool,
          requestId: 'call-fixture',
          status: 'succeeded',
          parts: [{ kind: 'text', text: 'result' }],
        },
      ],
    },
    declares(second),
    { role: 'assistant', parts: [{ kind: 'text', text: 'Ready' }] },
  ];
}

describe('recorded wire compatibility', () => {
  it.each(fixtures.fixtures)('sends the bytes the vendored adapter sent: $name', async (fixture) => {
    const { provider } = fixture;
    const recording = recordingFetch(responses[provider]('done'));
    const adapter = factories[provider]({ resolveApiKey: () => 'fixture-secret', fetch: recording.fetch });
    const parameters: ProjectInput['configuration']['parameters'] =
      fixture.mode === 'reasoning'
        ? {
            reasoning:
              provider === 'openai'
                ? { mode: 'off' }
                : { mode: 'budget', tokens: provider === 'anthropic' ? 2048 : 512 },
          }
        : {};
    const configuration = { target: { provider, model: models[provider] }, parameters };
    // The collected fixtures predate native loading on these models: the capabilities say so.
    const reported = await adapter.capabilitiesFor(configuration);
    const capabilities: ModelCapabilities = {
      ...reported,
      appendsTools: fixture.mode === 'append-only' && reported.appendsTools,
    };
    const messages = conversation();
    const tail: Message = { role: 'user', parts: [{ kind: 'text', text: 'Append é 🙂' }] };
    const input = { configuration, tools: [lookupTool, second], capabilities };
    for (const retained of [messages, [...messages, tail]]) {
      const request = await adapter.project({ ...input, messages: retained });
      expect(resultOf(await collectFrames(adapter, request)).finish).toBe('complete');
    }
    expect(recording.requests.map((request) => createHash('sha256').update(request.body).digest('hex'))).toEqual(
      fixture.wireDigests,
    );
  });
});
