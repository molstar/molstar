/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ModelCapabilities, ModelConfiguration, ModelParameters } from '../../messages/index.js';
import { type ModelLimits, modelInfo } from '../models.js';

/** A model's limits as the application gives them, over the catalog's. */
export type ModelLimitOverride = ModelLimits;
export type ModelLimitOverrides = Readonly<Record<string, ModelLimitOverride>>;

export const defaultMaxOutputTokens = 4096;

/** Whether the catalog (`../models.ts`) says a model appends tools where they load. */
export function appendsTools(provider: string, model: string): boolean {
  return modelInfo(provider, model)?.quirks.appendsTools === true;
}

/** A model's capabilities: the application's override or the catalog's limits (under `catalog`, the provider's name
 * there; false for an endpoint that is not a provider's own), plus what the provider can do. Custom and gateway models
 * need explicit limits. */
export function resolveProviderCapabilities(
  configuration: ModelConfiguration,
  provider: string,
  overrides: ModelLimitOverrides,
  features: Pick<ModelCapabilities, 'appendsTools' | 'cacheScope' | 'outputLimitIncludesReasoning'>,
  catalog: string | false = provider,
): ModelCapabilities {
  if (configuration.target.provider !== provider) {
    throw new Error(`${provider} adapter cannot project provider ${configuration.target.provider}.`);
  }
  const model = configuration.target.model;
  const selected = overrides[model] ?? (catalog ? modelInfo(catalog, model)?.limits : undefined);
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
