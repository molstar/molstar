/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createAnthropic } from '@ai-sdk/anthropic';
import type { Json, JsonObject, ModelParameters, Provider, ToolDefinition } from '../../messages/index.js';
import { UnrepresentableInput } from '../../messages/index.js';
import type { AttachmentRule } from '../common/attachments.js';
import {
  appendsTools,
  defaultMaxOutputTokens,
  type ModelLimitOverride,
  resolveProviderCapabilities,
} from '../common/capabilities.js';
import { estimateTool, type TokenEstimate } from '../common/estimate.js';
import {
  type ContextItem,
  type Converter,
  convertMessage,
  type ReasoningPart,
  type SdkMessage,
  systemText,
  uniqueTools,
  wireName,
} from '../common/messages.js';
import { createProvider, type Layout } from '../common/provider.js';
import { type ReplayState, type StreamPart, streamTextDelta } from '../common/stream.js';

export interface AnthropicProviderOptions {
  readonly resolveApiKey: () => string | Promise<string>;
  /** A calibrated text-token profile for this provider (default: four characters a token). */
  readonly tokenEstimate?: TokenEstimate;
  readonly fetch?: typeof globalThis.fetch;
  readonly headers?: Readonly<Record<string, string>>;
  /** Explicit limits for custom/gateway targets or application-specific restrictions. */
  readonly modelLimits?: Readonly<Record<string, ModelLimitOverride>>;
  /** Five-minute cache breakpoints on the system prompt, the top level and the last immediate tool; default true. */
  readonly cache?: boolean;
}

const images = new Set(['image/jpeg', 'image/png', 'image/gif', 'image/webp']);

const attachments: AttachmentRule = {
  name: 'anthropic',
  accepts: (type, placement) =>
    (placement === 'user' || placement === 'tool-result') && (images.has(type) || type === 'application/pdf'),
  textAsFile: (placement) => placement === 'user',
};

function parameters(parameters: ModelParameters): void {
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

const cacheControl = { type: 'ephemeral', ttl: '5m' } as const;

/** The appended layout (`tool_addition`): the tool list is the head's inventory, deferred except what loads before the
 * conversation starts; a later load becomes a `tool_addition` in a system message where it loads. */
function appended(c: Converter, items: readonly ContextItem[], inventory: readonly ToolDefinition[]): Layout {
  const initial = new Set<string>();
  const system: string[] = [];
  const messages: SdkMessage[] = [];
  let started = false;
  let inlineText: string[] = [];
  let changes: string[] = [];

  function flush(nextRole?: string): void {
    if (inlineText.length === 0 && changes.length === 0) return;
    const previous = messages.at(-1)?.role;
    if ((previous !== 'user' && previous !== 'tool') || (nextRole !== undefined && nextRole !== 'assistant')) {
      throw new UnrepresentableInput(
        'Anthropic inline updates require a complete user/tool-result boundary before an assistant or request end.',
      );
    }
    // The SDK otherwise hoists plain text into top-level system when no initial system exists.
    // A native tool change forces the legal inline path.
    if (inlineText.length > 0 && changes.length === 0 && system.length === 0) {
      throw new UnrepresentableInput('Anthropic in-place system text requires an initial system prompt with this SDK.');
    }
    messages.push({
      role: 'system',
      content: inlineText.join('\n\n'),
      ...(changes.length === 0
        ? {}
        : {
            providerOptions: {
              anthropic: { toolChanges: changes.map((toolName) => ({ type: 'tool_addition', toolName })) },
            },
          }),
    });
    inlineText = [];
    changes = [];
  }

  for (const item of items) {
    if (item.kind === 'load') {
      const id = item.definition.tool.id;
      if (!started) initial.add(id);
      else changes.push(wireName(id));
      continue;
    }
    if (item.message.role === 'system') {
      const text = systemText(c, item.message);
      if (text.length > 0) (started ? inlineText : system).push(text);
      continue;
    }
    // A tool load takes effect at the end of its user turn (later user messages in that turn go
    // first); in-place system text stays where it is written.
    if (item.message.role === 'assistant' || inlineText.length > 0) flush(item.message.role);
    started = true;
    messages.push(...convertMessage(c, item.message));
  }
  flush();
  return {
    system: system.join('\n\n'),
    messages,
    tools: uniqueTools(inventory).map((tool) => (initial.has(tool.id) ? tool : { ...tool, deferLoading: true })),
  };
}

const record = (value: Json | undefined): JsonObject | undefined =>
  value !== null && typeof value === 'object' && !Array.isArray(value) ? (value as JsonObject) : undefined;

/** Signed or redacted thinking kept on a call (`{ anthropic: { thinking: [...] } }`), as the reasoning before it. */
function thinking(options: Json | undefined): ReasoningPart[] {
  const blocks = record(record(options)?.anthropic)?.thinking;
  if (!Array.isArray(blocks)) return [];
  return blocks.flatMap((b): ReasoningPart[] => {
    const block = record(b);
    if (block === undefined) return [];
    if (typeof block.signature === 'string')
      return [
        {
          type: 'reasoning',
          text: typeof block.text === 'string' ? block.text : '',
          providerOptions: { anthropic: { signature: block.signature } },
        },
      ];
    if (typeof block.redactedData === 'string')
      return [{ type: 'reasoning', text: '', providerOptions: { anthropic: { redactedData: block.redactedData } } }];
    return [];
  });
}

interface ThinkingBlock {
  text: string;
  signature?: string;
  redactedData?: string;
}

/** Thinking streamed before a call goes onto that call, for its replay. */
function collect(): ReplayState {
  const blocks = new Map<string, ThinkingBlock>();
  let pending: ThinkingBlock[] = [];
  const keep = (part: StreamPart) => {
    const id = part.id ?? String(blocks.size);
    let block = blocks.get(id);
    if (!block) {
      block = { text: '' };
      blocks.set(id, block);
      pending.push(block);
    }
    if (part.type === 'reasoning-delta') block.text += streamTextDelta(part);
    const meta = part.providerMetadata?.anthropic;
    if (typeof meta?.signature === 'string' && meta.signature.length > 0) block.signature = meta.signature;
    if (typeof meta?.redactedData === 'string' && meta.redactedData.length > 0) block.redactedData = meta.redactedData;
  };
  return {
    reasoning(part) {
      if (part.type === 'reasoning-start' || part.type === 'reasoning-delta') keep(part);
    },
    toolCall(call) {
      const replay = pending.filter((block) => block.signature !== undefined || block.redactedData !== undefined);
      pending = [];
      if (replay.length === 0) return call;
      return {
        ...call,
        providerOptions: {
          ...call.providerOptions,
          anthropic: {
            thinking: replay.map((block) => ({
              text: block.text,
              ...(block.signature === undefined ? {} : { signature: block.signature }),
              ...(block.redactedData === undefined ? {} : { redactedData: block.redactedData }),
            })),
          },
        },
      };
    },
    finish: (calls) => calls,
  };
}

/** Anthropic over the AI SDK. Caches by prefix with breakpoints, so a cache scope id is ignored. */
export function createAnthropicProvider(options: AnthropicProviderOptions): Provider {
  const cache = options.cache !== false;
  return createProvider({
    provider: 'anthropic',
    label: 'Anthropic',
    revision: 'anthropic-ai-sdk-v4',
    capabilities: (configuration) =>
      resolveProviderCapabilities(configuration, 'anthropic', options.modelLimits ?? {}, {
        appendsTools: appendsTools('anthropic', configuration.target.model),
        cacheScope: false,
        outputLimitIncludesReasoning: true,
      }),
    resolveApiKey: options.resolveApiKey,
    parameters,
    attachments,
    replay: { reasoning: thinking, collect },
    native: { supports: (model) => appendsTools('anthropic', model), project: appended },
    call({ system, reasoning, tools }) {
      const lastImmediate = tools.findLastIndex((tool) => tool.deferLoading !== true);
      const anthropic: JsonObject = {
        ...(cache ? { cacheControl } : {}),
        ...(reasoning === undefined
          ? {}
          : reasoning.mode === 'budget'
            ? { thinking: { type: 'enabled', budgetTokens: reasoning.tokens } }
            : { thinking: { type: 'adaptive' }, ...(reasoning.effort ? { effort: reasoning.effort } : {}) }),
      };
      return {
        ...(system.length === 0
          ? {}
          : {
              system: cache
                ? { role: 'system', content: system, providerOptions: { anthropic: { cacheControl } } }
                : system,
            }),
        ...(Object.keys(anthropic).length === 0 ? {} : { providerOptions: { anthropic } }),
        toolOptions(tool, index) {
          const breakpoint = cache && index === lastImmediate;
          if (tool.deferLoading !== true && !breakpoint) return undefined;
          return {
            anthropic: {
              ...(tool.deferLoading === true ? { deferLoading: true } : {}),
              ...(breakpoint ? { cacheControl } : {}),
            },
          };
        },
      };
    },
    estimate: {
      // Claude injects tool-use instructions: its published one-tool count_tokens example
      // records 403 tokens, beyond schema text alone. This fixed allowance is a heuristic.
      overhead: (layout) => (layout.tools.length > 0 ? 313 : 0),
      disclosure(message, tools) {
        if (message.role !== 'system') return 0;
        const changes = record(message.providerOptions?.anthropic as Json | undefined)?.toolChanges;
        return (Array.isArray(changes) ? changes : []).reduce<number>((sum, change) => {
          const tool = tools.get(String(record(change)?.toolName));
          if (tool === undefined) throw new Error('Projected disclosure references an undeclared tool.');
          return sum + estimateTool(tool, options.tokenEstimate);
        }, 0);
      },
      ...(options.tokenEstimate ? { profile: options.tokenEstimate } : {}),
    },
    languageModel: (apiKey, payload) =>
      createAnthropic({
        apiKey,
        ...(options.fetch === undefined ? {} : { fetch: options.fetch }),
        ...(options.headers === undefined ? {} : { headers: { ...options.headers } }),
      })(payload.model),
  });
}
