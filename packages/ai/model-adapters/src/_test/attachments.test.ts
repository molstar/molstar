/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Message, ModelVisiblePart } from '@molstar/ai-agent';
import { UnrepresentableInput } from '@molstar/ai-agent';
import { bytesDigest, inlineAttachment } from '../common/attachments.js';
import { createAnthropicModelAdapter, createGoogleModelAdapter, createOpenAIModelAdapter } from '../index.js';
import {
  anthropicSse,
  collectFrames,
  googleSse,
  lookupTool,
  payloadOf,
  project,
  projectInput,
  recordingFetch,
} from './support.js';

function file(mediaType: string, content = 'fixture'): ModelVisiblePart {
  return {
    kind: 'reference',
    reference: { uri: `data:${mediaType};base64,${btoa(content)}`, mediaType, name: 'fixture.file' },
  };
}

it.each([
  { provider: 'anthropic', model: 'claude-sonnet-5-5', create: createAnthropicModelAdapter },
  { provider: 'google', model: 'gemini-3.8-flash', create: createGoogleModelAdapter },
  { provider: 'openai', model: 'gpt-6-luna', create: createOpenAIModelAdapter },
] as const)(
  'classifies unsupported $provider content separately from corrupt retained bytes',
  async ({ provider, model, create }) => {
    const adapter = create({ resolveApiKey: () => 'fixture' });
    for (const message of [
      { role: 'user', parts: [file('application/octet-stream')] },
      { role: 'system', parts: [file('image/png')] },
      {
        role: 'user',
        parts: [
          { kind: 'tool-result', requestId: 'c', tool: lookupTool.tool, status: 'failed', parts: [file('image/png')] },
        ],
      },
    ] satisfies readonly Message[]) {
      await expect(project(adapter, projectInput(provider, model, { messages: [message] }))).rejects.toBeInstanceOf(
        UnrepresentableInput,
      );
    }
    const corrupt: Message = {
      role: 'user',
      parts: [{ kind: 'reference', reference: { uri: 'data:image/png;base64,YQ==', digest: 'sha256:bad' } }],
    };
    await expect(project(adapter, projectInput(provider, model, { messages: [corrupt] }))).rejects.not.toBeInstanceOf(
      UnrepresentableInput,
    );
  },
);

describe('inline attachment receiving boundary', () => {
  it('recognizes arbitrary media and normalizes binary percent data without provider assumptions', async () => {
    expect(await inlineAttachment({ uri: 'data:application/x-fixture,%00%FF%2B' })).toMatchObject({
      mediaType: 'application/x-fixture',
      data: 'AP8r',
      byteLength: 3,
      digest: await bytesDigest(new Uint8Array([0, 255, 43])),
    });
    expect(await bytesDigest(new TextEncoder().encode('abc'))).toBe(
      'sha256:ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad',
    );
    expect(
      await inlineAttachment({ uri: 'https://example.test/file.pdf', mediaType: 'application/pdf' }),
    ).toBeUndefined();
    expect(await inlineAttachment({ uri: 'data:,hello%20world' })).toMatchObject({
      mediaType: 'text/plain',
      byteLength: 11,
    });
  });

  it.each([
    [{ uri: 'data:image/png;base64' }, /comma/],
    [{ uri: 'data:image/png;base64,%%%' }, /base64/],
    [{ uri: 'data:image/png;base64,Zh==' }, /base64/],
    [{ uri: 'data:image/png,%GG' }, /percent/],
    [{ uri: 'data:image/png;charset=latin1,abc' }, /encoding/],
    [{ uri: 'data:image/png;base64,YQ==', mediaType: 'application/pdf' }, /media type/],
    [{ uri: 'data:text/plain;charset=utf-8,abc', mediaType: 'text/plain;charset=us-ascii' }, /charset/],
    [{ uri: 'data:image/png;base64,YQ==', byteLength: 2 }, /byteLength/],
    [{ uri: 'data:image/png;base64,YQ==', digest: 'sha256:bad' }, /digest/],
  ] as const)('rejects invalid retained inline bytes %#', async (reference, error) => {
    await expect(inlineAttachment(reference)).rejects.toThrow(error);
  });
});

for (const provider of ['anthropic', 'google'] as const) {
  const model = provider === 'anthropic' ? 'claude-sonnet-4-5' : 'gemini-3-flash-preview';
  const create = provider === 'anthropic' ? createAnthropicModelAdapter : createGoogleModelAdapter;
  const response = provider === 'anthropic' ? anthropicSse('done') : googleSse('done');
  describe(`${provider} attachment SDK projection`, () => {
    it('retains files-only input and preserves ordered text, documents and ordinary references', async () => {
      const recorded = recordingFetch(response);
      const adapter = create({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
      const request = await project(
        adapter,
        projectInput(provider, model, {
          messages: [
            {
              role: 'user',
              parts: [
                { kind: 'text', text: 'before' },
                file('image/png'),
                { kind: 'text', text: 'between' },
                file('application/pdf'),
                file('text/plain', 'document text'),
                { kind: 'reference', reference: { uri: 'resource:immutable' } },
              ],
            },
            { role: 'user', parts: [file('image/jpeg')] },
          ],
        }),
      );
      await collectFrames(adapter, request);
      const body = JSON.parse(recorded.requests[0]?.body ?? '{}');
      if (provider === 'anthropic') {
        expect(body.messages[0].content.map((part: { type: string }) => part.type)).toEqual([
          'text',
          'image',
          'text',
          'document',
          'document',
          'text',
          'image',
        ]);
        expect(body.messages[0].content[4].source).toEqual({
          type: 'text',
          media_type: 'text/plain',
          data: 'document text',
        });
        expect(body.messages[0].content[1].source.data).toBe(btoa('fixture'));
      } else {
        expect(body.contents[0].parts).toEqual([
          { text: 'before' },
          { inlineData: { mimeType: 'image/png', data: btoa('fixture') } },
          { text: 'between' },
          { inlineData: { mimeType: 'application/pdf', data: btoa('fixture') } },
          { text: 'document text' },
          { text: '[ref resource:immutable]' },
        ]);
        expect(body.contents[1].parts).toEqual([{ inlineData: { mimeType: 'image/jpeg', data: btoa('fixture') } }]);
      }
    });

    it('keeps a mixed result group adjacent to its calls before visible text and files', async () => {
      const recorded = recordingFetch(response);
      const adapter = create({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
      const messages: Message[] = [
        { role: 'user', parts: [{ kind: 'text', text: 'start' }] },
        {
          role: 'assistant',
          parts: [{ kind: 'tool-call', requestId: 'toolu_fixture', tool: lookupTool.tool, arguments: { id: 'x' } }],
        },
        {
          role: 'user',
          parts: [
            { kind: 'text', text: 'follow up' },
            file('image/jpeg', 'outside'),
            {
              kind: 'tool-result',
              requestId: 'toolu_fixture',
              tool: lookupTool.tool,
              status: 'succeeded',
              parts: [
                { kind: 'text', text: 'result' },
                file('image/png', 'inside'),
                file('application/pdf', 'pdf'),
                file('text/markdown', '# result'),
              ],
            },
          ],
        },
      ];
      await collectFrames(adapter, await project(adapter, projectInput(provider, model, { messages })));
      const body = JSON.parse(recorded.requests[0]?.body ?? '{}');
      if (provider === 'anthropic') {
        const content = body.messages[2].content;
        expect(content.map((part: { type: string }) => part.type)).toEqual(['tool_result', 'text', 'image']);
        expect(content[0].content.map((part: { type: string }) => part.type)).toEqual([
          'text',
          'image',
          'document',
          'text',
        ]);
        expect(content[0].content[1].source.data).toBe(btoa('inside'));
      } else {
        const result = body.contents[2].parts[0].functionResponse;
        expect(result.id).toBe('toolu_fixture');
        expect(result.parts).toEqual([
          { inlineData: { mimeType: 'image/png', data: btoa('inside') } },
          { inlineData: { mimeType: 'application/pdf', data: btoa('pdf') } },
        ]);
        expect(body.contents[3].parts).toEqual([
          { text: 'follow up' },
          { inlineData: { mimeType: 'image/jpeg', data: btoa('outside') } },
        ]);
      }
    });

    it('rejects unsupported MIME/role combinations and invalid text before invocation', async () => {
      const recorded = recordingFetch(response);
      const adapter = create({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
      for (const [role, part] of [
        ['user', file('application/octet-stream')],
        ['system', file('image/png')],
        ['user', { kind: 'reference', reference: { uri: 'data:text/plain;base64,/w==' } }],
      ] as const) {
        await expect(
          project(adapter, projectInput(provider, model, { messages: [{ role, parts: [part] }] })),
        ).rejects.toThrow(/attachment/);
      }
      expect(recorded.requests).toHaveLength(0);
    });

    it('retains byte evidence of an attachment in the request', async () => {
      const recorded = recordingFetch(response);
      const adapter = create({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
      const request = await project(
        adapter,
        projectInput(provider, model, { messages: [{ role: 'user', parts: [file('image/png')] }] }),
      );
      expect(JSON.stringify(payloadOf(request).messages)).toContain(
        await bytesDigest(new TextEncoder().encode('fixture')),
      );
      await collectFrames(adapter, request);
      expect(recorded.requests).toHaveLength(1);
    });

    it('rejects binary failed-result attachments rather than losing their error status', async () => {
      const adapter = create({ resolveApiKey: () => 'fixture' });
      const failed: Message = {
        role: 'user',
        parts: [
          { kind: 'tool-result', requestId: 'c', tool: lookupTool.tool, status: 'failed', parts: [file('image/png')] },
        ],
      };
      const failure = await project(adapter, projectInput(provider, model, { messages: [failed] })).catch(
        (error: unknown) => error,
      );
      expect(failure).toBeInstanceOf(UnrepresentableInput);
      expect((failure as Error).message).toMatch(/failed tool-result/);
    });

    it('sizes decoded text and modality bytes rather than base64 prompt characters', async () => {
      const adapter = create({ resolveApiKey: () => 'fixture' });
      const inline = await project(
        adapter,
        projectInput(provider, model, {
          messages: [{ role: 'user', parts: [file('text/plain', 'text '.repeat(1000))] }],
        }),
      );
      const text = await project(
        adapter,
        projectInput(provider, model, {
          messages: [{ role: 'user', parts: [{ kind: 'text', text: 'text '.repeat(1000) }] }],
        }),
      );
      expect(Math.abs(inline.inputTokens - text.inputTokens)).toBeLessThan(128);
      const image = await project(
        adapter,
        projectInput(provider, model, { messages: [{ role: 'user', parts: [file('image/png', 'a'.repeat(12000))] }] }),
      );
      expect(image.inputTokens).toBeGreaterThan(4000);
      const percent = await project(
        adapter,
        projectInput(provider, model, {
          messages: [
            {
              role: 'user',
              parts: [{ kind: 'reference', reference: { uri: `data:image/png,${'%61'.repeat(12000)}` } }],
            },
          ],
        }),
      );
      expect(percent.inputTokens).toBe(image.inputTokens);
    });
  });
}

it('maps Google audio/video and assistant attachments and rejects Claude audio', async () => {
  const recorded = recordingFetch(googleSse('done'));
  const adapter = createGoogleModelAdapter({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
  await collectFrames(
    adapter,
    await project(
      adapter,
      projectInput('google', 'gemini-3-flash-preview', {
        messages: [
          { role: 'user', parts: [file('audio/wav'), file('video/mp4')] },
          { role: 'assistant', parts: [file('image/png')] },
        ],
      }),
    ),
  );
  const body = JSON.parse(recorded.requests[0]?.body ?? '{}');
  expect(body.contents[0].parts.map((part: { inlineData: { mimeType: string } }) => part.inlineData.mimeType)).toEqual([
    'audio/wav',
    'video/mp4',
  ]);
  expect(body.contents[1].parts[0].inlineData.mimeType).toBe('image/png');
  const claude = createAnthropicModelAdapter({ resolveApiKey: () => 'fixture' });
  await expect(
    project(
      claude,
      projectInput('anthropic', 'claude-sonnet-4-5', { messages: [{ role: 'user', parts: [file('audio/wav')] }] }),
    ),
  ).rejects.toThrow(/audio\/wav/);
  await expect(
    project(
      claude,
      projectInput('anthropic', 'claude-sonnet-4-5', { messages: [{ role: 'assistant', parts: [file('image/png')] }] }),
    ),
  ).rejects.toThrow(/assistant/);
});

it('rejects binary Google result attachments on targets whose SDK relocates them', async () => {
  const adapter = createGoogleModelAdapter({ resolveApiKey: () => 'fixture' });
  await expect(
    project(
      adapter,
      projectInput('google', 'gemini-2.5-flash', {
        messages: [
          {
            role: 'user',
            parts: [
              {
                kind: 'tool-result',
                requestId: 'c',
                tool: lookupTool.tool,
                status: 'succeeded',
                parts: [file('image/png')],
              },
            ],
          },
        ],
      }),
    ),
  ).rejects.toThrow(/tool-result/);
});
