/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type {
  ModelCapabilities,
  ModelConfiguration,
  ModelFrame,
  ProjectedRequest,
  ProjectInput,
  Provider,
} from '../../messages/index.js';
import { createOpenAIProvider } from '../../providers/openai/provider.js';
import { scriptedModel } from '../../testing/index.js';
import { AgentSession } from '../index.js';
import { guidance, type Toy, toy, toyToolkit } from './helpers.js';

/** A session over one provider lookup, the default model `target`. */
const sessionOver = (providerFor: (id: string) => Provider | undefined, target: { provider: string; model: string }) =>
  new AgentSession<Toy>({ models: { providerFor, default: target }, toolkit: toyToolkit(), api: toy(), guidance });

describe('the model’s limits from a real provider (models.md, "The providers")', () => {
  it('clamps the output it asks for to a model whose maximum is below the session’s default', async () => {
    // An own OpenAI-compatible endpoint with explicit limits: 8k output, under the session's 16k.
    const real = createOpenAIProvider({
      resolveApiKey: () => 'unused',
      baseURL: 'http://localhost:11434/v1',
      api: 'chat',
      provider: 'own-lab',
      modelLimits: { llama: { contextWindowTokens: 128_000, maxOutputTokens: 8192 } },
    });
    // The real provider's capabilities and projection; the scripted model's stream (no network).
    const scripted = scriptedModel([{ text: 'Hello.' }]);
    const projected: number[] = [];
    const lab: Provider = {
      capabilities: (configuration) => real.capabilities(configuration),
      project: async (input) => {
        const request = await real.project(input);
        projected.push((request.payload as { maxOutputTokens: number }).maxOutputTokens);
        return request;
      },
      invoke: (input) => scripted.provider.invoke(input),
    };
    const session = sessionOver((id) => (id === 'own-lab' ? lab : undefined), { provider: 'own-lab', model: 'llama' });
    const r = await session.turn({ prompt: 'Hi.' });
    expect(r.error).toBeUndefined();
    expect(r.status).toBe('completed');
    expect(projected).toEqual([8192]);
    expect(session.contextWindow()).toBe(128_000);
  });

  it('keeps the methods of a provider made by a class', async () => {
    const scripted = scriptedModel([{ text: 'From a class.' }]);
    class ClassProvider implements Provider {
      private readonly inner = scripted.provider;
      capabilities(configuration: ModelConfiguration): Promise<ModelCapabilities> {
        return this.inner.capabilities(configuration);
      }
      project(input: ProjectInput): Promise<ProjectedRequest> {
        return this.inner.project(input);
      }
      invoke(input: { readonly request: ProjectedRequest; readonly signal: AbortSignal }): AsyncIterable<ModelFrame> {
        return this.inner.invoke(input);
      }
    }
    const provider = new ClassProvider();
    const session = sessionOver(() => provider, { provider: 'scripted', model: 'scripted-1' });
    const r = await session.turn({ prompt: 'Hi.' });
    expect(r).toMatchObject({ status: 'completed', text: 'From a class.' });
  });
});

describe('a provider the host replaces between turns', () => {
  it('runs the next turn on the new provider (an own endpoint whose address changed)', async () => {
    const before = scriptedModel([{ text: 'From the old address.' }]);
    const after = scriptedModel([{ text: 'From the new address.' }]);
    let current = before.provider;
    const session = sessionOver(() => current, { provider: 'own-lab', model: 'llama' });
    expect((await session.turn({ prompt: 'One.' })).text).toBe('From the old address.');
    current = after.provider;
    const r = await session.turn({ prompt: 'Two.' });
    expect(r).toMatchObject({ status: 'completed', text: 'From the new address.' });
    expect([before.invocations.length, after.invocations.length]).toEqual([1, 1]);
  });
});
