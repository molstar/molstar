/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, test } from 'vitest';
import type { JsonObject, Message, ProjectInput, ToolDefinition } from '../../messages/index.js';
import { createAnthropicModelAdapter, createGoogleModelAdapter, createOpenAIModelAdapter } from '../index.js';
import {
  anthropicSse,
  collectFrames,
  declares,
  googleSse,
  lookupTool,
  openaiSse,
  payloadOf,
  project,
  projectInput,
  recordingFetch,
} from './support.js';

const user: Message = { role: 'user', parts: [{ kind: 'text', text: 'Choose the appropriate tool.' }] };
const deferred: ToolDefinition = {
  ...lookupTool,
  tool: { id: 'application.extra' },
  name: 'application.extra',
  description: 'A complete deferred definition. '.repeat(1_000),
};
const gateway = { contextWindowTokens: 1_000_000, maxOutputTokens: 4096 };
const limits = { 'gpt-gateway': gateway };

// Native cases load a tool where its part stands; collected cases take the tools from `tools`.
const cases = [
  { name: 'anthropic native', provider: 'anthropic', model: 'claude-sonnet-5-5', native: true },
  { name: 'anthropic collected', provider: 'anthropic', model: 'claude-sonnet-5', native: false },
  { name: 'google collected', provider: 'google', model: 'gemini-3.8-flash', native: false },
  { name: 'openai native', provider: 'openai', model: 'gpt-6-luna', native: true },
  { name: 'openai collected', provider: 'openai', model: 'gpt-gateway', native: false },
] as const;
const factories = {
  anthropic: createAnthropicModelAdapter,
  google: createGoogleModelAdapter,
  openai: createOpenAIModelAdapter,
};
const responses = { anthropic: anthropicSse, google: googleSse, openai: openaiSse };

function setup(test: (typeof cases)[number]) {
  const recorded = recordingFetch(responses[test.provider]('ok'));
  const adapter = factories[test.provider]({
    resolveApiKey: () => 'fixture',
    fetch: recorded.fetch,
    modelLimits: limits,
  });
  const input = (messages: Message[], tools: ToolDefinition[]): Omit<ProjectInput, 'capabilities'> =>
    projectInput(test.provider, test.model, { messages, tools });
  const body = (index: number): JsonObject => JSON.parse(recorded.requests[index]?.body ?? '{}');
  return { recorded, adapter, input, body };
}
const history: Message[] = [
  user,
  {
    role: 'assistant',
    parts: [{ kind: 'tool-call', requestId: 'historical-call', tool: deferred.tool, arguments: {} }],
  },
  {
    role: 'user',
    parts: [
      {
        kind: 'tool-result',
        requestId: 'historical-call',
        tool: deferred.tool,
        status: 'succeeded',
        parts: [{ kind: 'text', text: 'Historical result' }],
      },
    ],
  },
];

describe('deferred tool definitions', () => {
  it.each(cases.flatMap((test) => [true, false].map((declared) => ({ ...test, declared }))))(
    'retains completed history with no loaded tools: $name, declared=$declared',
    async (test) => {
      const { adapter, input, body } = setup(test);
      const declaration = test.declared ? [declares(deferred, true)] : [];
      const compacted = await project(adapter, input([...declaration, ...history], []));
      expect(payloadOf(compacted).messages).toEqual(
        expect.arrayContaining([
          expect.objectContaining({ toolCalls: [expect.objectContaining({ requestId: 'historical-call' })] }),
          expect.objectContaining({ toolResults: [expect.objectContaining({ requestId: 'historical-call' })] }),
        ]),
      );
      const reloaded = await project(adapter, input([...declaration, ...history, declares(deferred)], [deferred]));
      expect(reloaded.inputTokens).toBeGreaterThan(compacted.inputTokens + 1000);
      await collectFrames(adapter, compacted);
      await collectFrames(adapter, reloaded);
      if (test.provider === 'openai') {
        const wire = body(0).input as JsonObject[];
        expect(wire).toEqual(
          expect.arrayContaining([
            expect.objectContaining({ type: 'function_call', call_id: 'historical-call' }),
            expect.objectContaining({
              type: 'function_call_output',
              call_id: 'historical-call',
              output: 'Historical result',
            }),
          ]),
        );
        expect(body(0).tools).toBeUndefined();
        if (test.native) {
          expect(wire.some((item) => item.type === 'additional_tools')).toBe(false);
          expect((body(1).input as JsonObject[]).at(-1)).toMatchObject({ type: 'additional_tools' });
        }
      }
    },
  );

  it.each(cases)('projects $name replacement heads as declared-only, then loads again', async (test) => {
    const { adapter, input, body, recorded } = setup(test);
    const summary: Message = {
      role: 'system',
      parts: [{ kind: 'text', text: 'Prior tool work has been summarized.' }],
    };
    const declaration = declares(deferred, true);
    const compactedMessages = [summary, declaration, user];
    const loaded = await project(adapter, input([summary, declares(deferred), user], [deferred]));
    const compacted = await project(adapter, input(compactedMessages, []));
    const reloaded = await project(adapter, input([...compactedMessages, declares(deferred)], [deferred]));
    expect(loaded.inputTokens).toBeGreaterThan(compacted.inputTokens + 1000);
    expect(reloaded.inputTokens).toBeGreaterThan(compacted.inputTokens + 1000);
    for (const request of [loaded, compacted, reloaded]) await collectFrames(adapter, request);
    expect(recorded.requests).toHaveLength(3);
    if (test.provider === 'anthropic' && test.native) {
      expect((body(1).tools as JsonObject[])[0]).toMatchObject({ defer_loading: true });
      expect(JSON.stringify(body(1).messages)).not.toContain('tool_addition');
      expect(JSON.stringify(body(2).messages)).toContain('tool_addition');
    } else if (test.provider === 'openai' && test.native) {
      expect(JSON.stringify(body(1).input)).not.toContain('additional_tools');
      expect(JSON.stringify(body(2).input)).toContain('additional_tools');
    } else {
      const compactedTools = body(1).tools;
      expect(compactedTools === undefined || (Array.isArray(compactedTools) && compactedTools.length === 0)).toBe(true);
    }
  });

  it('keeps new deferred declarations out of the active earlier prompt and sizing', async () => {
    const test = cases[0];
    const { adapter, input, body } = setup(test);
    const head = [declares(lookupTool), user];
    const first = await project(adapter, input(head, [lookupTool]));
    const declared = await project(adapter, input([...head, declares(deferred, true)], [lookupTool]));
    expect(declared.inputTokens).toBe(first.inputTokens);
    expect(declared.digest).not.toBe(first.digest);
    await collectFrames(adapter, first);
    await collectFrames(adapter, declared);
    expect(body(1).messages).toEqual(body(0).messages);
    expect(body(1).system).toEqual(body(0).system);
    expect((body(1).tools as JsonObject[])[0]).toEqual((body(0).tools as JsonObject[])[0]);
    expect((body(1).tools as JsonObject[])[1]).toMatchObject({
      name: 'application_extra',
      defer_loading: true,
      description: deferred.description,
      input_schema: deferred.inputSchema,
    });
  });

  it('loads full schemas at their position and ignores identical repeat loads', async () => {
    const { adapter, input, body } = setup(cases[0]);
    const head = [declares(lookupTool), user];
    const declaredOnly = await project(adapter, input([...head, declares(deferred, true)], [lookupTool]));
    const loadedMessages = [...head, declares(deferred)];
    const loaded = await project(adapter, input(loadedMessages, [lookupTool, deferred]));
    const repeated = await project(adapter, input([...loadedMessages, declares(deferred)], [lookupTool, deferred]));
    expect(loaded.inputTokens).toBeGreaterThan(declaredOnly.inputTokens + 1_000);
    expect(repeated).toEqual(loaded);
    await collectFrames(adapter, loaded);
    expect((body(0).messages as JsonObject[]).at(-1)).toEqual({
      role: 'system',
      content: [{ type: 'tool_addition', tool: { type: 'tool_reference', name: 'application_extra' } }],
    });
  });

  it('loads a tool given without a loading part before the conversation starts', async () => {
    const { adapter, input } = setup(cases[0]);
    const request = await project(adapter, input([user], [lookupTool]));
    const tools = payloadOf(request).tools as JsonObject[];
    expect(tools).toHaveLength(1);
    expect(tools[0]?.deferLoading).toBeUndefined();
  });

  it.each([cases[0], cases[3]])(
    'loads what the head loads, whatever the request may call (a nested run keeps the prefix): $name',
    async (test) => {
      const { adapter, input } = setup(test);
      const head = [declares(lookupTool), user, declares(deferred)];
      const parent = payloadOf(await project(adapter, input(head, [lookupTool, deferred])));
      const nested = payloadOf(await project(adapter, input(head, [])));
      expect(nested.tools).toEqual(parent.tools);
      expect(nested.messages).toEqual(parent.messages);
      expect((nested.tools as JsonObject[]).some((tool) => tool.deferLoading !== true)).toBe(true);
    },
  );

  it('maps a call to a declared but unloaded tool to its exact declared id', async () => {
    const unloaded: ToolDefinition = { ...deferred, tool: { id: 'application.complex/id:extra' } };
    const recorded = recordingFetch(openaiSse('', { name: 'application_complex_id_extra', input: '{"id":"x"}' }));
    const adapter = createOpenAIModelAdapter({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
    const request = await project(
      adapter,
      projectInput('openai', 'gpt-6-luna', {
        messages: [declares(lookupTool), declares(unloaded, true), user],
        tools: [lookupTool],
      }),
    );
    const frames = await collectFrames(adapter, request);
    expect(frames.find((frame) => frame.kind === 'result')).toMatchObject({
      message: { parts: [{ kind: 'tool-call', tool: unloaded.tool }] },
    });
    const wire = recorded.requests[0]?.body ?? '';
    expect(wire).not.toContain(unloaded.description);
    expect(wire).not.toContain('application_complex_id_extra"');
  });

  it('replays the retained declaration and load bytes despite later factory caching preferences', async () => {
    const firstWire = recordingFetch(anthropicSse('ok'));
    const first = createAnthropicModelAdapter({ resolveApiKey: () => 'fixture', fetch: firstWire.fetch });
    const request = await project(
      first,
      projectInput('anthropic', 'claude-sonnet-5-5', {
        messages: [declares(lookupTool), user, declares(deferred)],
        tools: [lookupTool, deferred],
      }),
    );
    const retained = JSON.parse(JSON.stringify(request));
    const replayWire = recordingFetch(anthropicSse('ok'));
    const replay = createAnthropicModelAdapter({
      resolveApiKey: () => 'fixture',
      fetch: replayWire.fetch,
      cache: false,
    });
    await collectFrames(first, request);
    await collectFrames(replay, retained);
    expect(JSON.parse(replayWire.requests[0]?.body ?? '{}')).toEqual(JSON.parse(firstWire.requests[0]?.body ?? '{}'));
  });
});
