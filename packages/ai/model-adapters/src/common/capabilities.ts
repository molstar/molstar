/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelCapabilities, ModelConfiguration, ModelParameters } from '@molstar/ai-agent';

export interface ModelLimitOverride {
  readonly contextWindowTokens: number;
  readonly maxInputTokens?: number;
  readonly maxOutputTokens: number;
  readonly defaultMaxOutputTokens?: number;
}

export type ModelLimits = Readonly<Record<string, ModelLimitOverride>>;

export const defaultMaxOutputTokens = 4096;

/** A model's capabilities: the bundled limits or the application's override, plus what the
 * provider can do. Custom and gateway models need explicit limits. */
export function resolveProviderCapabilities(
  configuration: ModelConfiguration,
  provider: string,
  limits: ModelLimits,
  overrides: ModelLimits,
  features: Pick<ModelCapabilities, 'appendsTools' | 'cacheScope' | 'outputLimitIncludesReasoning'>,
): ModelCapabilities {
  if (configuration.target.provider !== provider) {
    throw new Error(`${provider} adapter cannot project provider ${configuration.target.provider}.`);
  }
  const model = configuration.target.model;
  const selected = overrides[model] ?? limits[model];
  if (selected === undefined) {
    throw new Error(`Unknown ${provider} model ${model}; supply a modelLimits override with explicit token limits.`);
  }
  const capabilities: ModelCapabilities = {
    ...selected,
    ...features,
    defaultMaxOutputTokens:
      selected.defaultMaxOutputTokens ?? Math.min(defaultMaxOutputTokens, selected.maxOutputTokens),
  };
  validateCapabilities(capabilities, configuration.parameters);
  return capabilities;
}

/** Whether a limit snapshot is usable and holds the requested output. */
export function validateCapabilities(capabilities: ModelCapabilities, parameters: ModelParameters): void {
  for (const [name, value] of Object.entries({
    contextWindowTokens: capabilities.contextWindowTokens,
    maxInputTokens: capabilities.maxInputTokens,
    maxOutputTokens: capabilities.maxOutputTokens,
    defaultMaxOutputTokens: capabilities.defaultMaxOutputTokens,
    requestedMaxOutputTokens: parameters.maxOutputTokens,
  })) {
    if (value !== undefined && (!Number.isSafeInteger(value) || value <= 0)) {
      throw new Error(`${name} must be a positive safe integer token count.`);
    }
  }
  if (
    capabilities.defaultMaxOutputTokens !== undefined &&
    capabilities.maxOutputTokens !== undefined &&
    capabilities.defaultMaxOutputTokens > capabilities.maxOutputTokens
  ) {
    throw new Error('Default output limit exceeds model maximum.');
  }
  const output = parameters.maxOutputTokens ?? capabilities.defaultMaxOutputTokens;
  if (output !== undefined && capabilities.maxOutputTokens !== undefined && output > capabilities.maxOutputTokens) {
    throw new Error(`Requested output limit ${output} exceeds model maximum ${capabilities.maxOutputTokens}.`);
  }
}
