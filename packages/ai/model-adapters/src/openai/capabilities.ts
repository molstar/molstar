/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelCapabilities, ModelConfiguration, ModelParameters } from '@molstar/ai-agent';
import { appendsTools, type ModelLimitOverrides, resolveProviderCapabilities } from '../common/capabilities.js';

/** Positional full tool additions (`additional_tools`), where the catalog says the model has them. */
export function supportsOpenAIAppendOnly(model: string): boolean {
  return appendsTools('openai', model);
}

/** OpenAI maps a cache scope id to the Responses `promptCacheKey`. Another endpoint (`own` false)
 * knows none of OpenAI's models, and chat completions have neither tool additions nor a cache key. */
export function providerCapabilities(
  configuration: ModelConfiguration,
  overrides: ModelLimitOverrides = {},
  endpoint: { provider: string; own: boolean; api: 'responses' | 'chat' } = {
    provider: 'openai',
    own: true,
    api: 'responses',
  },
): ModelCapabilities {
  const responses = endpoint.api === 'responses';
  return resolveProviderCapabilities(
    configuration,
    endpoint.provider,
    overrides,
    {
      appendsTools: endpoint.own && responses && supportsOpenAIAppendOnly(configuration.target.model),
      cacheScope: responses,
      outputLimitIncludesReasoning: true,
    },
    endpoint.own ? 'openai' : false,
  );
}

export function assertSupportedParameters(parameters: ModelParameters, _model: string): void {
  if (parameters.temperature !== undefined || parameters.topP !== undefined) {
    throw new Error("OpenAI Responses temperature and topP are not supported by this adapter's reasoning targets.");
  }
  const reasoning = parameters.reasoning;
  if (reasoning?.mode === 'budget') {
    throw new Error('OpenAI Responses reasoning token budgets are not supported; use adaptive effort.');
  }
  if (
    reasoning?.mode === 'adaptive' &&
    reasoning.effort !== undefined &&
    !['none', 'low', 'medium', 'high', 'xhigh', 'max'].includes(reasoning.effort)
  ) {
    throw new Error('openai reasoning effort is unsupported.');
  }
}
