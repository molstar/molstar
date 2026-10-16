/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, ModelAdapter, ModelCapabilities, ModelConfiguration, ProjectedRequest } from '@molstar/ai-agent';
import { digest, ModelError } from '@molstar/ai-agent';
import type { streamText } from 'ai';
import {
  buildProjectedPayload,
  type MeasurementPolicy,
  measureProjected,
  type ProjectedPayload,
  type ProjectionPolicy,
} from './payload.js';
import { invokeProjectedRequest, modelError, type StreamPolicy } from './stream.js';

export interface ProviderSpec {
  /** The provider id requests are made for (a target's provider): a vendor, or an endpoint of the same dialect. */
  readonly provider: string;
  readonly label: string;
  /** The identity of the projection and replay rules, independent of package versions. */
  readonly revision: string;
  readonly projection: ProjectionPolicy;
  readonly measurement: MeasurementPolicy;
  readonly stream: StreamPolicy;
  readonly capabilities: (configuration: ModelConfiguration) => ModelCapabilities;
  readonly resolveApiKey: () => string | Promise<string>;
  readonly languageModel: (
    apiKey: string,
    model: string,
    request: ProjectedRequest,
  ) => Parameters<typeof streamText>[0]['model'];
}

/** A part's digest, short: enough to tell parts apart in one conversation. */
const short = (value: Json) => digest(value).slice(0, 12);

/** The payload's parts in the order providers cache them: tools, system, then each message. */
function segmentsOf(payload: ProjectedPayload): string[] {
  return [
    short(payload.tools as unknown as Json),
    short(payload.system),
    ...payload.messages.map((m) => short(m as Json)),
  ];
}

/** The adapter every provider shares: projection is a pure function of its input, invocation streams
 * the retained payload through the AI SDK. */
export function createProviderAdapter(spec: ProviderSpec): ModelAdapter {
  return {
    capabilitiesFor: async (configuration) => spec.capabilities(configuration),
    async project(input) {
      const payload = await buildProjectedPayload(input, spec.provider, spec.revision, spec.projection);
      const json = payload as unknown as Json;
      return {
        payload: json,
        digest: digest(json),
        adapterRevision: spec.revision,
        inputTokens: measureProjected(payload, spec.measurement),
        segments: segmentsOf(payload),
      };
    },
    async *invoke(input) {
      const apiKey = await spec.resolveApiKey();
      if (apiKey.length === 0) {
        throw new ModelError(`${spec.label} API key is required for live model inference.`, { effect: 'none' });
      }
      try {
        const model = (input.request.payload as { model?: unknown } | null)?.model;
        if (typeof model !== 'string' || model.length === 0) {
          throw new ModelError(`Projected ${spec.label} request is missing a model id.`, { effect: 'none' });
        }
        yield* invokeProjectedRequest({
          request: input.request,
          signal: input.signal,
          expectedRevision: spec.revision,
          expectedProvider: spec.provider,
          policy: spec.stream,
          languageModel: spec.languageModel(apiKey, model, input.request),
        });
      } catch (error) {
        throw modelError(error, apiKey);
      }
    },
  };
}
