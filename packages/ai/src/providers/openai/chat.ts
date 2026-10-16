/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, JsonObject, ModelParameters, ModelVisiblePart } from '../../messages/index.js';
import {
  type AttachmentPlacement,
  projectVisibleContent as projectCommonContent,
  unsupportedAttachment,
} from '../common/attachments.js';
import type { ProjectionPolicy } from '../common/payload.js';
import { type ReasoningPart, type StreamPolicy, streamTextDelta } from '../common/stream.js';

// Chat completions, as OpenAI-compatible servers speak it (vLLM, SGLang, LiteLLM, Ollama, ...):
// tools collected on every request (no positional additions), the system prompt first, reasoning
// read from `reasoning_content` and sent back with the tool calls it led to.

const images = new Set(['image/jpeg', 'image/png', 'image/gif', 'image/webp']);

/** Text attachments as text; images only in the person's own messages (tool messages are text). */
function projectVisibleContent(
  parts: readonly ModelVisiblePart[],
  model: string,
  placement: AttachmentPlacement,
): Promise<readonly JsonObject[]> {
  return projectCommonContent(parts, model, placement, (attachment, model, placement, text): JsonObject => {
    if (text !== undefined) return { type: 'text', text, attachment: { ...attachment } };
    if (placement !== 'user' || !images.has(attachment.mediaType)) {
      throw unsupportedAttachment('openai-chat', model, attachment.mediaType, placement);
    }
    return { type: 'file', ...attachment };
  });
}

function assertSupportedParameters(parameters: ModelParameters, _model: string): void {
  const reasoning = parameters.reasoning;
  if (reasoning?.mode === 'budget') {
    throw new Error('Chat completions take a reasoning effort, not a token budget.');
  }
  if (
    reasoning?.mode === 'adaptive' &&
    reasoning.effort !== undefined &&
    !['none', 'minimal', 'low', 'medium', 'high'].includes(reasoning.effort)
  ) {
    throw new Error('Chat completions reasoning effort is unsupported.');
  }
}

export const chatProjection: ProjectionPolicy = {
  visible: projectVisibleContent,
  assertParameters: assertSupportedParameters,
};

/** `name` is the provider's options key (the SDK passes its unknown fields into the request body);
 * `body` holds fields the server takes on every request (e.g. LiteLLM's `cache`). */
export function chatStreamPolicy(name: string, body: JsonObject = {}): StreamPolicy {
  return {
    options(payload) {
      const effort = payload.reasoning?.mode === 'adaptive' ? payload.reasoning.effort : undefined;
      const fields: JsonObject = { ...body, ...(effort === undefined ? {} : { reasoningEffort: effort }) };
      return {
        ...(payload.system.length > 0 ? { system: payload.system } : {}),
        ...(Object.keys(fields).length > 0 ? { providerOptions: { [name]: fields } } : {}),
      };
    },
    inlineMessage(message) {
      if (message.role === 'system') throw new Error('Projected provider does not support inline system messages.');
      return undefined;
    },
    reasoning: chatReasoning,
    createReplay() {
      // What the model reasoned before a batch of calls goes back with the batch's first call.
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
    },
  };
}

function chatReasoning(value: Json | undefined): ReasoningPart[] {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) return [];
  const chat = (value as JsonObject).chat;
  if (chat === null || typeof chat !== 'object' || Array.isArray(chat)) return [];
  const text = (chat as JsonObject).text;
  return typeof text === 'string' && text.length > 0 ? [{ type: 'reasoning', text, providerOptions: {} }] : [];
}
