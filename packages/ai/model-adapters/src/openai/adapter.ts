/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createOpenAI } from '@ai-sdk/openai';
import { createOpenAICompatible } from '@ai-sdk/openai-compatible';
import type { JsonObject, ModelAdapter } from '@molstar/ai-agent';
import { createProviderAdapter } from '../common/adapter.js';
import type { ModelLimitOverride } from '../common/capabilities.js';
import type { TokenEstimate } from '../common/estimate.js';
import { providerCapabilities } from './capabilities.js';
import { chatProjection, chatStreamPolicy } from './chat.js';
import { openAIDisclosureTransport } from './disclosure-transport.js';
import { measurementPolicy, projectionPolicy } from './payload.js';
import { streamPolicy } from './stream.js';

export interface OpenAIModelAdapterOptions {
  readonly resolveApiKey: () => string | Promise<string>;
  /** A calibrated text-token profile for this provider (default: four characters a token). */
  readonly tokenEstimate?: TokenEstimate;
  readonly fetch?: typeof globalThis.fetch;
  readonly modelLimits?: Readonly<Record<string, ModelLimitOverride>>;
  /** Another endpoint than OpenAI's own, OpenAI-compatible ("https://llm.ai.e-infra.cz/v1"). Its
   * models are known only by `modelLimits`. */
  readonly baseURL?: string;
  /** Which API it speaks: OpenAI's Responses (the default) or chat completions, which most
   * compatible servers speak natively. */
  readonly api?: 'responses' | 'chat';
  /** The provider requests are made for, as targets name it (default "openai"), and its label. */
  readonly provider?: string;
  readonly label?: string;
  /** Chat only: fields the server takes on every request (LiteLLM's `cache: { "no-cache": true }`). */
  readonly body?: JsonObject;
}

/** Stateless Responses projection; provider state is retained in the request/call replay evidence.
 * With `api: "chat"`, chat completions: tools collected, reasoning sent back with its calls. */
export function createOpenAIModelAdapter(options: OpenAIModelAdapterOptions): ModelAdapter {
  const provider = options.provider ?? 'openai';
  const api = options.api ?? 'responses';
  const own = options.baseURL === undefined;
  const fetch = options.fetch ?? globalThis.fetch;
  const common = {
    provider,
    label: options.label ?? 'OpenAI',
    capabilities: (configuration: Parameters<typeof providerCapabilities>[0]) =>
      providerCapabilities(configuration, options.modelLimits, { provider, own, api }),
    resolveApiKey: options.resolveApiKey,
  };
  const profile = options.tokenEstimate ? { profile: options.tokenEstimate } : {};
  if (api === 'chat') {
    return createProviderAdapter({
      ...common,
      revision: 'openai-chat-ai-sdk-v1',
      projection: chatProjection,
      measurement: profile,
      stream: chatStreamPolicy(provider, options.body),
      languageModel: (apiKey, model) =>
        createOpenAICompatible({
          name: provider,
          baseURL: options.baseURL ?? 'https://api.openai.com/v1',
          apiKey,
          fetch,
          includeUsage: true,
        }).chatModel(model),
    });
  }
  return createProviderAdapter({
    ...common,
    revision: 'openai-responses-ai-sdk-v3',
    projection: projectionPolicy,
    measurement: { ...measurementPolicy, ...profile },
    stream: streamPolicy,
    languageModel: (apiKey, model, request) =>
      createOpenAI({
        apiKey,
        ...(own ? {} : { baseURL: options.baseURL }),
        fetch: openAIDisclosureTransport(request, fetch),
      }).responses(model),
  });
}
