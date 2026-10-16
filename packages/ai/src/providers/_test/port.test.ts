/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, test } from 'vitest';
import type { ModelVisiblePart, ToolDefinition } from '../../messages/index.js';
import { type Json, type JsonObject, type Message, ModelError, UnrepresentableInput } from '../../messages/index.js';
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

const marker = 'DEFINITION-MARKER-7f3a';
const documented: ToolDefinition = {
  tool: { id: 'docs.search' },
  name: 'Docs search',
  description: `${marker} searches the documentation.`,
  inputSchema: { type: 'object', properties: { query: { type: 'string', description: marker } } },
};
const user: Message = { role: 'user', parts: [{ kind: 'text', text: 'Find it.' }] };
const guide: Message = {
  role: 'system',
  parts: [
    { kind: 'text', text: 'Be brief.' },
    { kind: 'tool-definition', definition: documented, deferred: false },
    { kind: 'tool-definition', definition: lookupTool, deferred: true },
  ],
};

const cases = [
  { name: 'anthropic native', provider: 'anthropic', model: 'claude-sonnet-5-5', create: createAnthropicModelAdapter },
  {
    name: 'anthropic collected',
    provider: 'anthropic',
    model: 'claude-haiku-4-5',
    create: createAnthropicModelAdapter,
  },
  { name: 'openai native', provider: 'openai', model: 'gpt-6-luna', create: createOpenAIModelAdapter },
  { name: 'google collected', provider: 'google', model: 'gemini-3.8-flash', create: createGoogleModelAdapter },
] as const;
const responses = { anthropic: anthropicSse, google: googleSse, openai: openaiSse };

/** Every string under a `text` or `system` key: the prompt text of a payload. */
function promptText(value: Json, out: string[] = []): string[] {
  if (Array.isArray(value)) for (const item of value) promptText(item, out);
  else if (value !== null && typeof value === 'object') {
    for (const [key, item] of Object.entries(value)) {
      if ((key === 'text' || key === 'system') && typeof item === 'string') out.push(item);
      else promptText(item, out);
    }
  }
  return out;
}

describe('tool-definition parts', () => {
  it.each(cases)('are never sent as text: $name', async ({ provider, model, create }) => {
    const recorded = recordingFetch(responses[provider]('ok'));
    const adapter = create({ resolveApiKey: () => 'fixture', fetch: recorded.fetch });
    const request = await project(
      adapter,
      projectInput(provider, model, { messages: [guide, user], tools: [documented] }),
    );
    const payload = payloadOf(request);
    expect(
      payload.system === 'Be brief.' || (payload.messages as JsonObject[]).some((m) => m.text === 'Be brief.'),
    ).toBe(true);
    for (const text of promptText(request.payload)) {
      expect(text).not.toContain(marker);
      expect(text).not.toContain('Docs search');
    }
    // The definition is declared to the provider, not said to the model.
    await collectFrames(adapter, request);
    expect(recorded.requests[0]?.body).toContain(marker);
    expect(recorded.requests[0]?.body).toContain('docs_search');
  });

  it.each(cases)('alone in a system message add no text: $name', async ({ provider, model, create }) => {
    const adapter = create({ resolveApiKey: () => 'fixture' });
    const withMessage = await project(
      adapter,
      projectInput(provider, model, { messages: [declares(documented), user], tools: [documented] }),
    );
    expect(payloadOf(withMessage).system ?? '').toBe('');
    expect(promptText(withMessage.payload).join('')).toBe('Find it.');
  });
});

describe('a tool loaded later', () => {
  const head = [guide, user];
  const later = [...head, declares(lookupTool)];

  it.each(cases.filter((test) => test.name.includes('native')))(
    'is appended where it loads, not hoisted: $name',
    async ({ provider, model, create }) => {
      const adapter = create({ resolveApiKey: () => 'fixture' });
      const first = payloadOf(
        await project(adapter, projectInput(provider, model, { messages: head, tools: [documented] })),
      );
      const second = payloadOf(
        await project(adapter, projectInput(provider, model, { messages: later, tools: [documented, lookupTool] })),
      );
      expect((second.messages as unknown[]).slice(0, (first.messages as unknown[]).length)).toEqual(first.messages);
      expect(second.messages as unknown[]).toHaveLength((first.messages as unknown[]).length + 1);
      const added = (second.messages as JsonObject[])[(first.messages as unknown[]).length];
      expect(JSON.stringify(added)).toMatch(provider === 'openai' ? /additional_tools/ : /tool_addition/);
      // The tools the SDK gets for the first request are a prefix of the second's.
      expect((second.tools as unknown[]).slice(0, (first.tools as unknown[]).length)).toEqual(first.tools);
    },
  );

  it.each(cases.filter((test) => test.name.includes('collected')))(
    'is collected from the available tools: $name',
    async ({ provider, model, create }) => {
      const adapter = create({ resolveApiKey: () => 'fixture' });
      const tools = async (messages: Message[], available: ToolDefinition[]) =>
        (
          payloadOf(await project(adapter, projectInput(provider, model, { messages, tools: available })))
            .tools as JsonObject[]
        ).map((tool) => tool.id);
      expect(await tools(head, [documented])).toEqual(['docs.search']);
      expect(await tools(later, [documented, lookupTool])).toEqual(['docs.search', 'fixture.lookup']);
      // Messages alone do not make a tool available.
      expect(await tools(later, [documented])).toEqual(['docs.search']);
    },
  );
});

describe('cacheScopeId', () => {
  const scoped = (provider: 'anthropic' | 'google' | 'openai', model: string, cacheScopeId?: string) => {
    const input = projectInput(provider, model);
    return {
      ...input,
      configuration: { ...input.configuration, ...(cacheScopeId === undefined ? {} : { cacheScopeId }) },
    };
  };

  it("changes OpenAI's payload and digest, and only OpenAI's", async () => {
    const openai = createOpenAIModelAdapter({ resolveApiKey: () => 'k' });
    const none = await project(openai, scoped('openai', 'gpt-6-luna'));
    const a = await project(openai, scoped('openai', 'gpt-6-luna', 'scope-a'));
    const b = await project(openai, scoped('openai', 'gpt-6-luna', 'scope-b'));
    expect(payloadOf(a).cacheScopeId).toBe('scope-a');
    expect(new Set([none.digest, a.digest, b.digest]).size).toBe(3);
    const anthropic = createAnthropicModelAdapter({ resolveApiKey: () => 'k' });
    const google = createGoogleModelAdapter({ resolveApiKey: () => 'k' });
    for (const [adapter, provider, model] of [
      [anthropic, 'anthropic', 'claude-sonnet-5-5'],
      [google, 'google', 'gemini-3.8-flash'],
    ] as const) {
      const plain = await project(adapter, scoped(provider, model));
      const withScope = await project(adapter, scoped(provider, model, 'scope-a'));
      expect(withScope.payload).toEqual(plain.payload);
      expect(withScope.digest).toBe(plain.digest);
    }
  });

  it('is ignored when the capabilities say the model has no cache scope', async () => {
    const openai = createOpenAIModelAdapter({ resolveApiKey: () => 'k' });
    const input = scoped('openai', 'gpt-6-luna', 'scope-a');
    const capabilities = { ...(await openai.capabilitiesFor(input.configuration)), cacheScope: false };
    expect(payloadOf(await openai.project({ ...input, capabilities })).cacheScopeId).toBeUndefined();
  });
});

describe('unrepresentable input', () => {
  const unrepresentable = async (run: Promise<unknown>) => {
    const failure = await run.then(
      () => undefined,
      (error: unknown) => error,
    );
    expect(failure).toBeInstanceOf(UnrepresentableInput);
    expect(failure).toMatchObject({ code: 'unrepresentable' });
  };
  const file = (mediaType: string): ModelVisiblePart => ({
    kind: 'reference',
    reference: { uri: `data:${mediaType};base64,YQ==`, mediaType },
  });

  it.each(cases)('is typed for an attachment the provider cannot take: $name', async ({ provider, model, create }) => {
    const adapter = create({ resolveApiKey: () => 'k' });
    await unrepresentable(
      project(
        adapter,
        projectInput(provider, model, { messages: [{ role: 'user', parts: [file('application/x-unknown')] }] }),
      ),
    );
    await unrepresentable(
      project(
        adapter,
        projectInput(provider, model, { messages: [{ role: 'system', parts: [file('image/png')] }, user] }),
      ),
    );
  });

  it('is typed for system text inside a conversation that collects it', async () => {
    const midway: Message = { role: 'system', parts: [{ kind: 'text', text: 'New rule.' }] };
    for (const test of cases.filter((candidate) => candidate.name.includes('collected'))) {
      const adapter = test.create({ resolveApiKey: () => 'k' });
      await unrepresentable(project(adapter, projectInput(test.provider, test.model, { messages: [user, midway] })));
    }
  });

  it('is typed for OpenAI guidance or a load while a call is unanswered', async () => {
    const adapter = createOpenAIModelAdapter({ resolveApiKey: () => 'k' });
    const call: Message = {
      role: 'assistant',
      parts: [{ kind: 'tool-call', requestId: 'c', tool: lookupTool.tool, arguments: {} }],
    };
    const guidance: Message = { role: 'system', parts: [{ kind: 'text', text: 'Late.' }] };
    await unrepresentable(project(adapter, projectInput('openai', 'gpt-6-luna', { messages: [user, call, guidance] })));
    await unrepresentable(
      project(
        adapter,
        projectInput('openai', 'gpt-6-luna', { messages: [user, call, declares(documented)], tools: [documented] }),
      ),
    );
  });

  it('leaves a bad parameter or a corrupt attachment a plain error', async () => {
    const adapter = createOpenAIModelAdapter({ resolveApiKey: () => 'k' });
    const input = projectInput('openai', 'gpt-6-luna');
    const bad = await project(adapter, {
      ...input,
      configuration: { ...input.configuration, parameters: { temperature: 0.2 } },
    }).catch((error: unknown) => error);
    expect(bad).toBeInstanceOf(Error);
    expect(bad).not.toBeInstanceOf(UnrepresentableInput);
    const corrupt = await project(adapter, {
      ...input,
      messages: [
        {
          role: 'user',
          parts: [{ kind: 'reference', reference: { uri: 'data:image/png;base64,YQ==', digest: 'sha256:bad' } }],
        },
      ],
    }).catch((error: unknown) => error);
    expect(corrupt).toBeInstanceOf(Error);
    expect(corrupt).not.toBeInstanceOf(UnrepresentableInput);
    expect(corrupt).not.toBeInstanceOf(ModelError);
  });
});

describe('pure projection', () => {
  it.each(cases)('gives the same payload and digest for the same input: $name', async ({ provider, model, create }) => {
    const input = projectInput(provider, model, { messages: [guide, user], tools: [documented] });
    const first = await project(create({ resolveApiKey: () => 'k' }), input);
    const second = await project(create({ resolveApiKey: () => 'other-key' }), structuredClone(input));
    expect(second).toEqual(first);
    expect(first.digest).toMatch(/^[a-f0-9]{32}$/);
    expect(Number.isInteger(first.inputTokens) && first.inputTokens > 0).toBe(true);
  });

  it('changes the digest when the configuration, a message or the tools change', async () => {
    const adapter = createAnthropicModelAdapter({ resolveApiKey: () => 'k' });
    const input = projectInput('anthropic', 'claude-sonnet-5-5', { messages: [user], tools: [lookupTool] });
    const base = (await project(adapter, input)).digest;
    const changed = [
      { ...input, tools: [documented] },
      { ...input, messages: [{ role: 'user', parts: [{ kind: 'text', text: 'Other.' }] } satisfies Message] },
      { ...input, configuration: { ...input.configuration, parameters: { maxOutputTokens: 64 } } },
    ];
    for (const next of changed) expect((await project(adapter, next)).digest).not.toBe(base);
  });
});
