/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createGoogleGenerativeAI } from '@ai-sdk/google';
import type { ModelAdapter } from '../../messages/index.js';
import { createProviderAdapter } from '../common/adapter.js';
import type { ModelLimitOverride } from '../common/capabilities.js';
import type { TokenEstimate } from '../common/estimate.js';
import type { MeasurementPolicy, ProjectionPolicy } from '../common/payload.js';
import { projectVisibleContent } from './attachments.js';
import { assertSupportedParameters, providerCapabilities } from './capabilities.js';
import { streamPolicy } from './stream.js';

export interface GoogleModelAdapterOptions {
  readonly resolveApiKey: () => string | Promise<string>;
  /** A calibrated text-token profile for this provider (default: four characters a token). */
  readonly tokenEstimate?: TokenEstimate;
  readonly fetch?: typeof globalThis.fetch;
  readonly modelLimits?: Readonly<Record<string, ModelLimitOverride>>;
}

const projection: ProjectionPolicy = { visible: projectVisibleContent, assertParameters: assertSupportedParameters };
const measurement: MeasurementPolicy = {};

export function createGoogleModelAdapter(options: GoogleModelAdapterOptions): ModelAdapter {
  return createProviderAdapter({
    provider: 'google',
    label: 'Google',
    revision: 'google-ai-sdk-v3',
    projection,
    measurement: { ...measurement, ...(options.tokenEstimate ? { profile: options.tokenEstimate } : {}) },
    stream: streamPolicy,
    capabilities: (configuration) => providerCapabilities(configuration, options.modelLimits),
    resolveApiKey: options.resolveApiKey,
    languageModel: (apiKey, model) =>
      createGoogleGenerativeAI({ apiKey, ...(options.fetch === undefined ? {} : { fetch: options.fetch }) })(model),
  });
}
