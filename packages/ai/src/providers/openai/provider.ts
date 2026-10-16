/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createOpenAI } from '@ai-sdk/openai';
import { createOpenAICompatible } from '@ai-sdk/openai-compatible';
import type { Json, JsonObject, ModelParameters, Provider, ToolCallPart } from '../../messages/index.js';
import type { AttachmentRule } from '../common/attachments.js';
import { appendsTools, type ModelLimitOverride, resolveProviderCapabilities } from '../common/capabilities.js';
import type { TokenEstimate } from '../common/estimate.js';
import type { ReasoningPart } from '../common/messages.js';
import { createProvider, type ProviderSpec } from '../common/provider.js';
import { type ReplayState, type StreamPart, streamTextDelta } from '../common/stream.js';
import { additionsTransport, additionTokens, appended } from './native.js';

export interface OpenAIProviderOptions {
  readonly resolveApiKey: () => string | Promise<string>;
  /** A calibrated text-token profile for this provider (default: four characters a token). */
  readonly tokenEstimate?: TokenEstimate;
  readonly fetch?: typeof globalThis.fetch;
  readonly modelLimits?: Readonly<Record<string, ModelLimitOverride>>;
  /** Another endpoint than OpenAI's own, OpenAI-compatible ("https://llm.ai.e-infra.cz/v1"). Its
   * models are known only by `modelLimits`. */
  readonly baseURL?: string;
  /** Which API it speaks: OpenAI's Responses (the default) or chat completions, which most
   * compatible servers speak natively. */
  readonly api?: 'responses' | 'chat';
  /** The provider requests are made for, as targets name it (default "openai"), and its label. */
  readonly provider?: string;
  readonly label?: string;
  /** Chat only: fields the server takes on every request (LiteLLM's `cache: { "no-cache": true }`). */
  readonly body?: JsonObject;
}

const images = new Set(['image/jpeg', 'image/png', 'image/gif', 'image/webp']);
// Document formats Responses takes, through its native input_file; a default file name keeps the extension for it.
const documents: Readonly<Record<string, string>> = {
  'application/pdf': 'pdf',
  'application/msword': 'doc',
  'application/vnd.openxmlformats-officedocument.wordprocessingml.document': 'docx',
  'application/rtf': 'rtf',
  'application/vnd.oasis.opendocument.text': 'odt',
  'application/vnd.ms-powerpoint': 'ppt',
  'application/vnd.openxmlformats-officedocument.presentationml.presentation': 'pptx',
  'application/vnd.ms-excel': 'xls',
  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet': 'xlsx',
};

const responsesAttachments: AttachmentRule = {
  name: 'openai',
  accepts: (type, placement) =>
    (placement === 'user' || placement === 'tool-result') && (images.has(type) || documents[type] !== undefined),
  filename: (type) => (documents[type] === undefined ? undefined : `attachment.${documents[type]}`),
};

/** Chat completions: images only in the person's own messages (tool messages are text). */
const chatAttachments: AttachmentRule = {
  name: 'openai-chat',
  accepts: (type, placement) => placement === 'user' && images.has(type),
};

const efforts = (parameters: ModelParameters, known: readonly string[], api: string) => {
  const reasoning = parameters.reasoning;
  if (reasoning?.mode === 'budget')
    throw new Error(`${api} reasoning token budgets are not supported: it takes an effort, not a token budget.`);
  if (reasoning?.mode === 'adaptive' && reasoning.effort !== undefined && !known.includes(reasoning.effort)) {
    throw new Error(`${api} reasoning effort is unsupported.`);
  }
};

function responsesParameters(parameters: ModelParameters): void {
  if (parameters.temperature !== undefined || parameters.topP !== undefined) {
    throw new Error("OpenAI Responses temperature and topP are not supported by this provider's reasoning targets.");
  }
  efforts(parameters, ['none', 'low', 'medium', 'high', 'xhigh', 'max'], 'OpenAI Responses');
}

const record = (value: Json | undefined): JsonObject | undefined =>
  value !== null && typeof value === 'object' && !Array.isArray(value) ? (value as JsonObject) : undefined;

/** Encrypted reasoning kept on a call (`{ openai: { reasoning: [...] } }`), as the reasoning items before it. */
function responsesReasoning(options: Json | undefined): ReasoningPart[] {
  const blocks = record(record(options)?.openai)?.reasoning;
  if (!Array.isArray(blocks)) return [];
  return blocks.flatMap((b): ReasoningPart[] => {
    const block = record(b);
    if (typeof block?.itemId !== 'string' || typeof block.encryptedContent !== 'string') return [];
    return [
      {
        type: 'reasoning',
        text: typeof block.text === 'string' ? block.text : '',
        providerOptions: { openai: { itemId: block.itemId, reasoningEncryptedContent: block.encryptedContent } },
      },
    ];
  });
}

/** The call's own item id. */
function responsesCallOptions(options: Json | undefined) {
  const itemId = record(record(options)?.openai)?.itemId;
  return typeof itemId === 'string' ? { openai: { itemId } } : undefined;
}

interface ReasoningBlock {
  itemId: string;
  text: string;
  encryptedContent?: string;
}

/** Reasoning items go onto the call they precede; end events can bring the encrypted content after the summary. */
function responsesCollect(): ReplayState {
  const blocks = new Map<string, ReasoningBlock>();
  let pending: ReasoningBlock[] = [];
  const byCall = new Map<number, ReasoningBlock[]>();
  const keep = (part: StreamPart) => {
    const meta = part.providerMetadata?.openai;
    const itemId = typeof meta?.itemId === 'string' ? meta.itemId : part.id?.split(':')[0];
    if (!itemId) return;
    let block = blocks.get(itemId);
    if (!block) {
      block = { itemId, text: '' };
      blocks.set(itemId, block);
      pending.push(block);
    }
    if (part.type === 'reasoning-delta') block.text += streamTextDelta(part);
    if (typeof meta?.reasoningEncryptedContent === 'string') block.encryptedContent = meta.reasoningEncryptedContent;
  };
  return {
    reasoning: keep,
    toolCall(call, part, index) {
      if (pending.length > 0) byCall.set(index, pending);
      pending = [];
      const itemId = part.providerMetadata?.openai?.itemId;
      return typeof itemId === 'string' && itemId.length > 0
        ? { ...call, providerOptions: { openai: { itemId } } }
        : call;
    },
    finish(calls): readonly ToolCallPart[] {
      return calls.map((call, index) => {
        const kept = byCall.get(index)?.filter((block) => block.encryptedContent !== undefined);
        if (!kept?.length) return call;
        return {
          ...call,
          providerOptions: {
            ...call.providerOptions,
            openai: {
              ...record(call.providerOptions?.openai),
              reasoning: kept.map((block) => ({
                itemId: block.itemId,
                text: block.text,
                encryptedContent: block.encryptedContent ?? '',
              })),
            },
          },
        };
      });
    },
  };
}

/** What the model reasoned before a batch of calls goes back with the batch's first call (`{ chat: { text } }`). */
function chatCollect(): ReplayState {
  let pending = '';
  return {
    reasoning(part) {
      if (part.type === 'reasoning-delta') pending += streamTextDelta(part);
    },
    toolCall(call) {
      const text = pending;
      pending = '';
      return text.length === 0 ? call : { ...call, providerOptions: { ...call.providerOptions, chat: { text } } };
    },
    finish: (calls) => calls,
  };
}

function chatReasoning(options: Json | undefined): ReasoningPart[] {
  const text = record(record(options)?.chat)?.text;
  return typeof text === 'string' && text.length > 0 ? [{ type: 'reasoning', text, providerOptions: {} }] : [];
}

/** OpenAI over the AI SDK: stateless Responses (the replay is on the calls), or with `api: "chat"` chat completions at
 * any compatible endpoint (vLLM, SGLang, LiteLLM, Ollama, ...): tools collected, reasoning sent back with its calls.
 * Another endpoint (`baseURL`) knows none of OpenAI's models and never appends tools. */
export function createOpenAIProvider(options: OpenAIProviderOptions): Provider {
  const provider = options.provider ?? 'openai';
  const api = options.api ?? 'responses';
  const own = options.baseURL === undefined;
  const fetch = options.fetch ?? globalThis.fetch;
  const responses = api === 'responses';
  const common: Pick<ProviderSpec, 'provider' | 'label' | 'capabilities' | 'resolveApiKey'> = {
    provider,
    label: options.label ?? 'OpenAI',
    capabilities: (configuration) =>
      resolveProviderCapabilities(
        configuration,
        provider,
        options.modelLimits ?? {},
        {
          appendsTools: own && responses && appendsTools('openai', configuration.target.model),
          // Responses maps a cache scope id to `promptCacheKey`; chat completions have none.
          cacheScope: responses,
          outputLimitIncludesReasoning: true,
        },
        own ? 'openai' : false,
      ),
    resolveApiKey: options.resolveApiKey,
  };
  if (!responses) {
    const body = options.body ?? {};
    return createProvider({
      ...common,
      revision: 'openai-chat-ai-sdk-v2',
      parameters: (parameters) => efforts(parameters, ['none', 'minimal', 'low', 'medium', 'high'], 'Chat completions'),
      attachments: chatAttachments,
      replay: { reasoning: chatReasoning, collect: chatCollect },
      call({ system, reasoning }) {
        // `provider` is the options key: the SDK passes its unknown fields into the request body.
        const effort = reasoning?.mode === 'adaptive' ? reasoning.effort : undefined;
        const fields: JsonObject = { ...body, ...(effort === undefined ? {} : { reasoningEffort: effort }) };
        return {
          ...(system.length > 0 ? { system } : {}),
          ...(Object.keys(fields).length > 0 ? { providerOptions: { [provider]: fields } } : {}),
        };
      },
      ...(options.tokenEstimate ? { estimate: { profile: options.tokenEstimate } } : {}),
      languageModel: (apiKey, payload) =>
        createOpenAICompatible({
          name: provider,
          baseURL: options.baseURL ?? 'https://api.openai.com/v1',
          apiKey,
          fetch,
          includeUsage: true,
        }).chatModel(payload.model),
    });
  }
  return createProvider({
    ...common,
    revision: 'openai-responses-ai-sdk-v4',
    parameters: responsesParameters,
    // "Off" is effort none: otherwise the model reasons by its own default.
    reasoning: (parameters) =>
      parameters.reasoning?.mode === 'off' ? { mode: 'adaptive', effort: 'none' } : parameters.reasoning,
    attachments: responsesAttachments,
    replay: { reasoning: responsesReasoning, callOptions: responsesCallOptions, collect: responsesCollect },
    native: { supports: (model) => appendsTools('openai', model), project: appended },
    call: ({ system, reasoning, cacheScopeId }) => ({
      ...(system.length > 0 ? { system } : {}),
      providerOptions: {
        openai: {
          store: false,
          passThroughUnsupportedFiles: true,
          ...(cacheScopeId === undefined ? {} : { promptCacheKey: cacheScopeId }),
          ...(reasoning?.mode === 'adaptive' && reasoning.effort !== undefined
            ? { reasoningEffort: reasoning.effort }
            : {}),
        },
      },
    }),
    estimate: {
      // The appended layout's tool list serves parsing only; the schemas count at their additions.
      countsTools: (isAppended) => !isAppended,
      disclosure: additionTokens,
      ...(options.tokenEstimate ? { profile: options.tokenEstimate } : {}),
    },
    languageModel: (apiKey, payload) =>
      createOpenAI({
        apiKey,
        ...(own ? {} : { baseURL: options.baseURL }),
        fetch: additionsTransport(payload, fetch),
      }).responses(payload.model),
  });
}
