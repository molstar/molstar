/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelCapabilities, ModelConfiguration, ModelParameters } from '@molstar/ai-agent';
import {
  defaultMaxOutputTokens,
  type ModelLimitOverride,
  type ModelLimits,
  resolveProviderCapabilities,
} from '../common/capabilities.js';

// Verified exercised targets; custom/gateway models require explicit limits.
const limits: ModelLimits = {
  'claude-haiku-4-5': { contextWindowTokens: 200_000, maxOutputTokens: 64_000 },
  'claude-haiku-4-5-20251001': { contextWindowTokens: 200_000, maxOutputTokens: 64_000 },
  'claude-haiku-5-5': { contextWindowTokens: 1_000_000, maxOutputTokens: 128_000 },
  'claude-sonnet-4-5': { contextWindowTokens: 200_000, maxOutputTokens: 64_000 },
  'claude-sonnet-4-5-20250929': { contextWindowTokens: 200_000, maxOutputTokens: 64_000 },
  'claude-sonnet-5': { contextWindowTokens: 1_000_000, maxOutputTokens: 128_000 },
  'claude-sonnet-5-5': { contextWindowTokens: 1_000_000, maxOutputTokens: 128_000 },
};

/** Native mid-conversation tool changes (`tool_addition`). */
export function supportsAnthropicAppendOnly(model: string): boolean {
  return model === 'claude-sonnet-5-5';
}

/** Anthropic caches by prefix with breakpoints, so a cache scope id is ignored. */
export function providerCapabilities(
  configuration: ModelConfiguration,
  overrides: Readonly<Record<string, ModelLimitOverride>> = {},
): ModelCapabilities {
  return resolveProviderCapabilities(configuration, 'anthropic', limits, overrides, {
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
