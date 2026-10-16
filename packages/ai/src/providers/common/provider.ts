/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  Json,
  JsonObject,
  ModelCapabilities,
  ModelConfiguration,
  ModelParameters,
  Provider,
  ReasoningParameters,
  ToolDefinition,
} from '../../messages/index.js';
import { digest, ModelError } from '../../messages/index.js';
import type { streamText } from 'ai';
import type { AttachmentRule } from './attachments.js';
import { defaultMaxOutputTokens, validateCapabilities } from './capabilities.js';
import { estimateContent, estimateTextTokens, estimateTool, type TokenEstimate } from './estimate.js';
import {
  type ContextItem,
  type Converter,
  collected,
  orderContext,
  type ProviderOptions,
  type ReplayCodec,
  type SdkMessage,
  type WireTool,
} from './messages.js';
import { modelError, type ReplayState, streamPayload } from './stream.js';

// The provider every vendor shares (providers.md §2): `project` builds `streamText`'s own arguments as JSON, the
// retained payload, digested; `invoke` checks it is this provider's and streams it as it is.

export interface PayloadTool {
  readonly description: string;
  readonly inputSchema: JsonObject;
  readonly providerOptions?: ProviderOptions;
}

/** The retained request: `streamText`'s arguments (but the model, a factory's), what `invoke` checks, and the tool
 * ids that calls map back to. JSON, never the API key. */
export interface Payload {
  readonly provider: string;
  readonly revision: string;
  readonly model: string;
  /** Wire name → tool id. */
  readonly toolIds: Readonly<Record<string, string>>;
  /** Tools load where they load (the `appendsTools` layout). */
  readonly appended?: true;
  readonly system?:
    | string
    | { readonly role: 'system'; readonly content: string; readonly providerOptions?: ProviderOptions };
  readonly messages: readonly SdkMessage[];
  readonly allowSystemInMessages?: true;
  readonly tools?: Readonly<Record<string, PayloadTool>>;
  readonly providerOptions?: ProviderOptions;
  readonly maxOutputTokens: number;
  readonly temperature?: number;
  readonly topP?: number;
}

/** A layout of the conversation: the system prompt, the messages, the tool list. */
export interface Layout {
  readonly system: string;
  readonly messages: readonly SdkMessage[];
  readonly tools: readonly WireTool[];
}

export type Reasoning = Exclude<ReasoningParameters, { readonly mode: 'off' }>;

export interface ProviderSpec {
  /** The provider id requests are made for (a target's provider): a vendor, or an endpoint of the same dialect. */
  readonly provider: string;
  readonly label: string;
  /** The identity of the projection and replay rules, independent of package versions; raise it when they change. */
  readonly revision: string;
  readonly capabilities: (configuration: ModelConfiguration) => ModelCapabilities;
  readonly resolveApiKey: () => string | Promise<string>;
  /** Refuses (a plain error) parameters the model cannot take. */
  readonly parameters: (parameters: ModelParameters, model: string) => void;
  /** The reasoning a request sends (OpenAI sends "off" as effort none); default: as given. */
  readonly reasoning?: (parameters: ModelParameters) => ReasoningParameters | undefined;
  readonly attachments: AttachmentRule;
  readonly replay: ReplayCodec & { collect(): ReplayState };
  /** The appended layout (tools loaded where they load), for the models whose catalog entry says `appendsTools`. */
  readonly native?: {
    readonly supports: (model: string) => boolean;
    readonly project: (c: Converter, items: readonly ContextItem[], inventory: readonly ToolDefinition[]) => Layout;
  };
  /** The call's own fields: the system prompt's form, provider options (reasoning, caching), each tool's options. */
  readonly call: (input: {
    readonly system: string;
    readonly reasoning?: Reasoning;
    readonly cacheScopeId?: string;
    readonly tools: readonly WireTool[];
  }) => {
    readonly system?: Payload['system'];
    readonly providerOptions?: ProviderOptions;
    readonly toolOptions?: (tool: WireTool, index: number) => ProviderOptions | undefined;
  };
  readonly estimate?: {
    /** Tokens a provider adds for the request (Claude's tool-use instructions). */
    readonly overhead?: (layout: Layout) => number;
    /** Whether the tool list is prompt (false: OpenAI's appended layout sends it only to parse calls). */
    readonly countsTools?: (appended: boolean) => boolean;
    /** A tool loaded mid-conversation, where it loads. */
    readonly disclosure?: (message: SdkMessage, tools: ReadonlyMap<string, WireTool>) => number;
    readonly profile?: TokenEstimate;
  };
  readonly languageModel: (apiKey: string, payload: Payload) => Parameters<typeof streamText>[0]['model'];
}

/** A part's digest, short: enough to tell parts apart in one conversation. */
const short = (value: unknown) => digest(value as Json).slice(0, 12);

/** The input estimate for admission: prompt content, not transport JSON or output controls (approximate,
 * deterministic). */
function estimate(spec: ProviderSpec, layout: Layout, appended: boolean): number {
  const policy = spec.estimate ?? {};
  const text = (value: string) => estimateTextTokens(value, policy.profile);
  const byName = new Map(layout.tools.map((tool) => [tool.name, tool]));
  let tokens = text(layout.system) + (policy.overhead?.(layout) ?? 0);
  for (const message of layout.messages) {
    tokens += 8 + estimateContent(message.content, text) + (policy.disclosure?.(message, byName) ?? 0);
  }
  if (policy.countsTools?.(appended) ?? true) {
    for (const tool of layout.tools) if (tool.deferLoading !== true) tokens += estimateTool(tool, policy.profile);
  }
  return Math.max(1, tokens);
}

/** A failure before any request: the request is refused or malformed, and nothing was generated. */
const refused = (message: string) => new ModelError(message, { transient: false, effect: 'none' });

export function createProvider(spec: ProviderSpec): Provider {
  return {
    capabilities: async (configuration) => spec.capabilities(configuration),
    async project({ messages, tools, configuration, capabilities }) {
      const model = configuration.target.model;
      if (configuration.target.provider !== spec.provider) {
        throw new Error(`The ${spec.provider} provider cannot project provider ${configuration.target.provider}.`);
      }
      const { parameters } = configuration;
      validateCapabilities(capabilities, parameters);
      spec.parameters(parameters, model);
      const reasoning = spec.reasoning?.(parameters) ?? parameters.reasoning;
      const c = { model, attachments: spec.attachments, replay: spec.replay, provider: spec.provider };
      let layout: Layout;
      if (capabilities.appendsTools) {
        if (spec.native?.supports(model) !== true) {
          throw new Error(`${spec.provider} model ${model} does not support appended tool definitions.`);
        }
        const { items, inventory } = orderContext(messages, tools);
        layout = spec.native.project(c, items, inventory);
      } else layout = collected(c, messages, tools);
      const call = spec.call({
        system: layout.system,
        ...(reasoning !== undefined && reasoning.mode !== 'off' ? { reasoning } : {}),
        ...(capabilities.cacheScope && configuration.cacheScopeId !== undefined
          ? { cacheScopeId: configuration.cacheScopeId }
          : {}),
        tools: layout.tools,
      });
      const payload: Payload = {
        provider: spec.provider,
        revision: spec.revision,
        model,
        toolIds: Object.fromEntries(layout.tools.map((tool) => [tool.name, tool.id])),
        ...(capabilities.appendsTools ? { appended: true } : {}),
        ...(call.system === undefined ? {} : { system: call.system }),
        messages: layout.messages,
        ...(layout.messages.some((message) => message.role === 'system') ? { allowSystemInMessages: true } : {}),
        ...(layout.tools.length === 0
          ? {}
          : {
              tools: Object.fromEntries(
                layout.tools.map((tool, index) => {
                  const options = call.toolOptions?.(tool, index);
                  return [
                    tool.name,
                    {
                      description: tool.description,
                      inputSchema: tool.inputSchema,
                      ...(options === undefined ? {} : { providerOptions: options }),
                    },
                  ];
                }),
              ),
            }),
        ...(call.providerOptions === undefined ? {} : { providerOptions: call.providerOptions }),
        maxOutputTokens: parameters.maxOutputTokens ?? capabilities.defaultMaxOutputTokens ?? defaultMaxOutputTokens,
        ...(parameters.temperature !== undefined ? { temperature: parameters.temperature } : {}),
        ...(parameters.topP !== undefined ? { topP: parameters.topP } : {}),
      };
      const json = payload as unknown as Json;
      return {
        payload: json,
        digest: digest(json),
        revision: spec.revision,
        inputTokens: estimate(spec, layout, capabilities.appendsTools),
        // The parts in the order providers cache them: tools, system, then each message (harness.md §11).
        prefix: [short(payload.tools ?? {}), short(payload.system ?? ''), ...payload.messages.map(short)],
      };
    },
    async *invoke({ request, signal }) {
      const apiKey = await spec.resolveApiKey();
      if (apiKey.length === 0) {
        throw new ModelError(`${spec.label} API key is required for live model inference.`, { effect: 'none' });
      }
      try {
        if (request.revision !== spec.revision) throw refused('Projected request belongs to an incompatible revision.');
        if (request.digest !== digest(request.payload))
          throw refused('Projected request payload digest does not match the payload.');
        const payload = request.payload as unknown as Payload | null;
        if (
          payload === null ||
          typeof payload !== 'object' ||
          !Array.isArray(payload.messages) ||
          typeof payload.model !== 'string' ||
          payload.model.length === 0
        ) {
          throw refused(`Projected ${spec.label} request is malformed.`);
        }
        if (payload.provider !== spec.provider || payload.revision !== spec.revision) {
          throw refused('Projected request provider does not match this provider.');
        }
        yield* streamPayload({
          payload,
          signal,
          languageModel: spec.languageModel(apiKey, payload),
          replay: spec.replay.collect(),
        });
      } catch (error) {
        throw modelError(error, apiKey);
      }
    },
  };
}
