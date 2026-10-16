/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createGoogleGenerativeAI } from '@ai-sdk/google';
import type { Json, JsonObject, ModelParameters, Provider } from '../../messages/index.js';
import type { AttachmentRule } from '../common/attachments.js';
import { type ModelLimitOverride, resolveProviderCapabilities } from '../common/capabilities.js';
import type { TokenEstimate } from '../common/estimate.js';
import { createProvider } from '../common/provider.js';

export interface GoogleProviderOptions {
  readonly resolveApiKey: () => string | Promise<string>;
  /** A calibrated text-token profile for this provider (default: four characters a token). */
  readonly tokenEstimate?: TokenEstimate;
  readonly fetch?: typeof globalThis.fetch;
  readonly modelLimits?: Readonly<Record<string, ModelLimitOverride>>;
}

const media = new Set([
  'application/pdf',
  'image/png',
  'image/jpeg',
  'image/webp',
  'image/heic',
  'image/heif',
  'audio/wav',
  'audio/mp3',
  'audio/mpeg',
  'audio/aiff',
  'audio/aac',
  'audio/ogg',
  'audio/flac',
  'video/mp4',
  'video/mpeg',
  'video/quicktime',
  'video/avi',
  'video/x-flv',
  'video/mpg',
  'video/webm',
  'video/wmv',
  'video/3gpp',
]);

/** Gemini 3 and later (`gemini-3`, `gemini-3.8-flash`, ...): reasoning parameters, media in tool results. */
export const isGemini3 = (model: string): boolean => /^gemini-(?:[3-9]|[1-9]\d)(?:[.-]|$)/.test(model);

/** Media anywhere but the system prompt; in a tool result only from Gemini 3 (the SDK moves them elsewhere before). */
const attachments: AttachmentRule = {
  name: 'google',
  accepts: (type, placement, model) =>
    placement !== 'system' && media.has(type) && (placement !== 'tool-result' || isGemini3(model)),
};

function parameters(parameters: ModelParameters, model: string): void {
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

/** A call's thought signature (`{ google: { thoughtSignature } }`), sent back on the call. */
function callOptions(value: Json | undefined) {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) return undefined;
  const google = (value as Record<string, Json>).google;
  if (google === null || typeof google !== 'object' || Array.isArray(google)) return undefined;
  const thoughtSignature = (google as Record<string, Json>).thoughtSignature;
  return typeof thoughtSignature === 'string' && thoughtSignature.length > 0
    ? { google: { thoughtSignature } }
    : undefined;
}

/** Gemini over the AI SDK: tools collected, no cache scope (its `cachedContent` has a lifecycle of its own). */
export function createGoogleProvider(options: GoogleProviderOptions): Provider {
  return createProvider({
    provider: 'google',
    label: 'Google',
    revision: 'google-ai-sdk-v4',
    capabilities: (configuration) =>
      resolveProviderCapabilities(configuration, 'google', options.modelLimits ?? {}, {
        appendsTools: false,
        cacheScope: false,
        outputLimitIncludesReasoning: false,
      }),
    resolveApiKey: options.resolveApiKey,
    parameters,
    attachments,
    replay: {
      callOptions,
      collect: () => ({
        reasoning() {},
        toolCall(call, part) {
          const thoughtSignature = part.providerMetadata?.google?.thoughtSignature;
          return typeof thoughtSignature === 'string' && thoughtSignature.length > 0
            ? { ...call, providerOptions: { google: { thoughtSignature } } }
            : call;
        },
        finish: (calls) => calls,
      }),
    },
    call({ system, reasoning }) {
      const thinkingConfig: JsonObject | undefined =
        reasoning === undefined
          ? undefined
          : reasoning.mode === 'budget'
            ? { thinkingBudget: reasoning.tokens }
            : reasoning.effort
              ? { thinkingLevel: reasoning.effort }
              : {};
      return {
        ...(system.length > 0 ? { system } : {}),
        ...(thinkingConfig === undefined ? {} : { providerOptions: { google: { thinkingConfig } } }),
      };
    },
    ...(options.tokenEstimate ? { estimate: { profile: options.tokenEstimate } } : {}),
    languageModel: (apiKey, payload) =>
      createGoogleGenerativeAI({ apiKey, ...(options.fetch === undefined ? {} : { fetch: options.fetch }) })(
        payload.model,
      ),
  });
}
