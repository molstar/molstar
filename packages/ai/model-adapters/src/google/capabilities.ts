/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelCapabilities, ModelConfiguration, ModelParameters } from '@molstar/ai-agent';
import { type ModelLimitOverride, type ModelLimits, resolveProviderCapabilities } from '../common/capabilities.js';
import { isGemini3 } from './attachments.js';

const window = { contextWindowTokens: 1_114_112, maxInputTokens: 1_048_576, maxOutputTokens: 65_536 };
// Verified exercised targets; custom/gateway models require explicit limits.
const limits: ModelLimits = {
  'gemini-2.5-flash': window,
  'gemini-3-flash-preview': window,
  'gemini-3.8-flash': window,
};

/** Gemini collects its tools and has no cache scope mapping (its `cachedContent` has a lifecycle). */
export function providerCapabilities(
  configuration: ModelConfiguration,
  overrides: Readonly<Record<string, ModelLimitOverride>> = {},
): ModelCapabilities {
  return resolveProviderCapabilities(configuration, 'google', limits, overrides, {
    appendsTools: false,
    cacheScope: false,
    outputLimitIncludesReasoning: false,
  });
}

export function assertSupportedParameters(parameters: ModelParameters, model: string): void {
  const reasoning = parameters.reasoning;
  if (reasoning === undefined || reasoning.mode === 'off') return;
  if (!isGemini3(model)) throw new Error('Google reasoning requires a Gemini 3 or later model.');
  if (reasoning.mode === 'budget') {
    if (!Number.isSafeInteger(reasoning.tokens) || reasoning.tokens < 1)
      throw new Error('google reasoning budget is invalid.');
  } else if (reasoning.effort !== undefined && !['minimal', 'low', 'medium', 'high'].includes(reasoning.effort)) {
    throw new Error('google reasoning effort is unsupported.');
  }
}
