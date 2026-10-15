/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelCapabilities, ModelConfiguration, ModelParameters } from '@molstar/ai-agent';
import {
  appendsTools,
  defaultMaxOutputTokens,
  type ModelLimitOverrides,
  resolveProviderCapabilities,
} from '../common/capabilities.js';

/** Native mid-conversation tool changes (`tool_addition`), where the catalog says the model has them. */
export function supportsAnthropicAppendOnly(model: string): boolean {
  return appendsTools('anthropic', model);
}

/** Anthropic caches by prefix with breakpoints, so a cache scope id is ignored. */
export function providerCapabilities(
  configuration: ModelConfiguration,
  overrides: ModelLimitOverrides = {},
): ModelCapabilities {
  return resolveProviderCapabilities(configuration, 'anthropic', overrides, {
    appendsTools: supportsAnthropicAppendOnly(configuration.target.model),
    cacheScope: false,
    outputLimitIncludesReasoning: true,
  });
}

export function assertSupportedParameters(parameters: ModelParameters, _model: string): void {
  const reasoning = parameters.reasoning;
  if (reasoning === undefined || reasoning.mode === 'off') return;
  if (reasoning.mode === 'budget') {
    if (!Number.isSafeInteger(reasoning.tokens) || reasoning.tokens < 1024)
      throw new Error('anthropic reasoning budget is invalid.');
    if (reasoning.tokens >= (parameters.maxOutputTokens ?? defaultMaxOutputTokens))
      throw new Error('Anthropic reasoning budget must be below maxOutputTokens.');
  } else if (reasoning.effort !== undefined && !['low', 'medium', 'high', 'xhigh', 'max'].includes(reasoning.effort)) {
    throw new Error('anthropic reasoning effort is unsupported.');
  }
}
