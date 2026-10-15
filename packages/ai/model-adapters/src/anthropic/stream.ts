/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, JsonObject } from '@molstar/ai-agent';
import { type ReasoningPart, type StreamPart, type StreamPolicy, streamTextDelta } from '../common/stream.js';

export const streamPolicy: StreamPolicy = {
  options(payload) {
    return {
      ...(payload.system.length > 0
        ? {
            system:
              payload.cache === true
                ? {
                    role: 'system' as const,
                    content: payload.system,
                    providerOptions: { anthropic: { cacheControl: cacheControl() } },
                  }
                : payload.system,
          }
        : {}),
      ...(payload.reasoning !== undefined || payload.cache === true
        ? {
            providerOptions: {
              anthropic: {
                ...(payload.cache === true ? { cacheControl: cacheControl() } : {}),
                ...(payload.reasoning === undefined
                  ? {}
                  : payload.reasoning.mode === 'budget'
                    ? { thinking: { type: 'enabled', budgetTokens: payload.reasoning.tokens } }
                    : {
                        thinking: { type: 'adaptive' },
                        ...(payload.reasoning.effort ? { effort: payload.reasoning.effort } : {}),
                      }),
              },
            },
          }
        : {}),
    };
  },
  toolOptions(tool, index, payload) {
    const lastImmediate = payload.tools.findLastIndex((item) => item.deferLoading !== true);
    if (tool.deferLoading !== true && !(payload.cache === true && index === lastImmediate)) return undefined;
    return {
      anthropic: {
        ...(tool.deferLoading === true ? { deferLoading: true } : {}),
        ...(payload.cache === true && index === lastImmediate ? { cacheControl: cacheControl() } : {}),
      },
    };
  },
  inlineMessage(message, _index, _digest) {
    if (message.role === 'system')
      return {
        role: 'system',
        content: typeof message.text === 'string' ? message.text : '',
        ...(Array.isArray(message.toolChanges)
          ? {
              providerOptions: {
                anthropic: {
                  toolChanges: message.toolChanges.map((change) => ({
                    type: 'tool_addition',
                    toolName: String((change as JsonObject).toolName),
                  })),
                },
              },
            }
          : {}),
      };
    return undefined;
  },
  reasoning: anthropicThinking,
  createReplay() {
    const thinking = new Map<string, ThinkingBlock>();
    let pending: ThinkingBlock[] = [];
    return {
      reasoning(part) {
        if (part.type === 'reasoning-start' || part.type === 'reasoning-delta')
          keepThinking(thinking, part, (block) => pending.push(block));
      },
      toolCall(call) {
        const replay = pending.filter((block) => block.signature !== undefined || block.redactedData !== undefined);
        pending = [];
        return replay.length === 0
          ? call
          : {
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
  },
};

interface ThinkingBlock {
  text: string;
  signature?: string;
  redactedData?: string;
}

function keepThinking(blocks: Map<string, ThinkingBlock>, part: StreamPart, added: (b: ThinkingBlock) => void): void {
  const id = part.id ?? String(blocks.size);
  let block = blocks.get(id);
  if (!block) {
    block = { text: '' };
    blocks.set(id, block);
    added(block);
  }
  if (part.type === 'reasoning-delta') block.text += streamTextDelta(part);
  const meta = part.providerMetadata?.anthropic;
  if (typeof meta?.signature === 'string' && meta.signature.length > 0) block.signature = meta.signature;
  if (typeof meta?.redactedData === 'string' && meta.redactedData.length > 0) block.redactedData = meta.redactedData;
}

function anthropicThinking(value: Json | undefined): ReasoningPart[] {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) return [];
  const anthropic = (value as JsonObject).anthropic;
  if (anthropic === null || typeof anthropic !== 'object' || Array.isArray(anthropic)) return [];
  const blocks = (anthropic as JsonObject).thinking;
  if (!Array.isArray(blocks)) return [];
  return blocks.flatMap((b): ReasoningPart[] => {
    if (b === null || typeof b !== 'object' || Array.isArray(b)) return [];
    const block = b as JsonObject;
    const text = typeof block.text === 'string' ? block.text : '';
    if (typeof block.signature === 'string')
      return [{ type: 'reasoning' as const, text, providerOptions: { anthropic: { signature: block.signature } } }];
    if (typeof block.redactedData === 'string')
      return [
        { type: 'reasoning' as const, text: '', providerOptions: { anthropic: { redactedData: block.redactedData } } },
      ];
    return [];
  });
}

function cacheControl(): { type: 'ephemeral'; ttl: '5m' } {
  return { type: 'ephemeral', ttl: '5m' };
}
