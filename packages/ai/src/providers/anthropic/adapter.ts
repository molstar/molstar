/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createAnthropic } from '@ai-sdk/anthropic';
import type { ModelAdapter } from '../../messages/index.js';
import { createProviderAdapter } from '../common/adapter.js';
import type { ModelLimitOverride } from '../common/capabilities.js';
import type { TokenEstimate } from '../common/estimate.js';
import { providerCapabilities } from './capabilities.js';
import { measurementPolicy, projectionPolicy } from './payload.js';
import { streamPolicy } from './stream.js';

export interface AnthropicModelAdapterOptions {
  readonly resolveApiKey: () => string | Promise<string>;
  /** A calibrated text-token profile for this provider (default: four characters a token). */
  readonly tokenEstimate?: TokenEstimate;
  readonly fetch?: typeof globalThis.fetch;
  readonly headers?: Readonly<Record<string, string>>;
  /** Explicit limits for custom/gateway targets or application-specific restrictions. */
  readonly modelLimits?: Readonly<Record<string, ModelLimitOverride>>;
  /** Native five-minute fixed-prefix and moving caching controls; enabled by default. */
  readonly cache?: boolean;
}

export function createAnthropicModelAdapter(options: AnthropicModelAdapterOptions): ModelAdapter {
  return createProviderAdapter({
    provider: 'anthropic',
    label: 'Anthropic',
    revision: 'anthropic-ai-sdk-v3',
    projection: projectionPolicy(options.cache !== false),
    measurement: { ...measurementPolicy, ...(options.tokenEstimate ? { profile: options.tokenEstimate } : {}) },
    stream: streamPolicy,
    capabilities: (configuration) => providerCapabilities(configuration, options.modelLimits),
    resolveApiKey: options.resolveApiKey,
    languageModel: (apiKey, model) =>
      createAnthropic({
        apiKey,
        ...(options.fetch === undefined ? {} : { fetch: options.fetch }),
        ...(options.headers === undefined ? {} : { headers: { ...options.headers } }),
      })(model),
  });
}
