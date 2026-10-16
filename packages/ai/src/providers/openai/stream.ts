/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, JsonObject, ToolCallPart } from '../../messages/index.js';
import { type ReasoningPart, type StreamPart, type StreamPolicy, streamTextDelta } from '../common/stream.js';
import { openAIDisclosureMarker } from './disclosure-transport.js';

export const streamPolicy: StreamPolicy = {
  options(payload) {
    return {
      ...(payload.system.length > 0 ? { system: payload.system } : {}),
      providerOptions: {
        openai: {
          store: false,
          passThroughUnsupportedFiles: true,
          ...(payload.cacheScopeId === undefined ? {} : { promptCacheKey: payload.cacheScopeId }),
          ...(payload.reasoning?.mode === 'adaptive' && payload.reasoning.effort !== undefined
            ? { reasoningEffort: payload.reasoning.effort }
            : {}),
        },
      },
    };
  },
  inlineMessage(message, index, digest) {
    if (message.type === 'additional_tools') {
      if (message.role !== 'developer' || !Array.isArray(message.tools))
        throw new Error('Invalid projected OpenAI native tool disclosure.');
      return { role: 'system', content: openAIDisclosureMarker(digest, index) };
    }
    if (message.role === 'system')
      return { role: 'system', content: typeof message.text === 'string' ? message.text : '' };
    return undefined;
  },
  reasoning: openAIThinking,
  replayOptions(value) {
    const options = openAIReplayOptions(value);
    return options === undefined ? undefined : { openai: { ...options } };
  },
  createReplay() {
    const blocks = new Map<string, OpenAIReasoningBlock>();
    let pending: OpenAIReasoningBlock[] = [];
    const callReasoning = new Map<number, OpenAIReasoningBlock[]>();
    return {
      reasoning(part) {
        keepOpenAIReasoning(blocks, part, (block) => pending.push(block));
      },
      toolCall(call, part, index) {
        if (pending.length > 0) callReasoning.set(index, pending);
        pending = [];
        const itemId = part.providerMetadata?.openai?.itemId;
        return typeof itemId === 'string' && itemId.length > 0
          ? { ...call, providerOptions: { openai: { itemId } } }
          : call;
      },
      finish(calls): readonly ToolCallPart[] {
        // End events can supply encrypted content after the summary was streamed.
        return calls.map((call, index) => {
          const blocks = callReasoning.get(index)?.filter((block) => block.encryptedContent !== undefined);
          if (!blocks?.length) return call;
          const existing = call.providerOptions?.openai;
          return {
            ...call,
            providerOptions: {
              ...call.providerOptions,
              openai: {
                ...(existing !== null && typeof existing === 'object' && !Array.isArray(existing) ? existing : {}),
                reasoning: blocks.map((block) => ({
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
  },
};

interface OpenAIReasoningBlock {
  itemId: string;
  text: string;
  encryptedContent?: string;
}

function keepOpenAIReasoning(
  blocks: Map<string, OpenAIReasoningBlock>,
  part: StreamPart,
  added: (block: OpenAIReasoningBlock) => void,
): void {
  const meta = part.providerMetadata?.openai;
  const itemId = typeof meta?.itemId === 'string' ? meta.itemId : part.id?.split(':')[0];
  if (!itemId) return;
  let block = blocks.get(itemId);
  if (!block) {
    block = { itemId, text: '' };
    blocks.set(itemId, block);
    added(block);
  }
  if (part.type === 'reasoning-delta') block.text += streamTextDelta(part);
  if (typeof meta?.reasoningEncryptedContent === 'string') block.encryptedContent = meta.reasoningEncryptedContent;
}

function openAIOptions(value: Json | undefined): JsonObject | undefined {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) return undefined;
  const openai = (value as JsonObject).openai;
  return openai !== null && typeof openai === 'object' && !Array.isArray(openai) ? (openai as JsonObject) : undefined;
}

function openAIThinking(value: Json | undefined): ReasoningPart[] {
  const blocks = openAIOptions(value)?.reasoning;
  if (!Array.isArray(blocks)) return [];
  return blocks.flatMap((entry): ReasoningPart[] => {
    if (entry === null || typeof entry !== 'object' || Array.isArray(entry)) return [];
    if (typeof entry.itemId !== 'string' || typeof entry.encryptedContent !== 'string') return [];
    return [
      {
        type: 'reasoning',
        text: typeof entry.text === 'string' ? entry.text : '',
        providerOptions: {
          openai: { itemId: entry.itemId, reasoningEncryptedContent: entry.encryptedContent },
        },
      },
    ];
  });
}

function openAIReplayOptions(value: Json | undefined): { itemId: string } | undefined {
  const options = openAIOptions(value);
  return typeof options?.itemId === 'string' ? { itemId: options.itemId } : undefined;
}
