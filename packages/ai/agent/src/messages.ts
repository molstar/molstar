/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, JsonObject, Reference } from './values.js';

/** What a model sees of content: text, or a reference the adapter resolves. */
export type ModelVisiblePart =
  | { readonly kind: 'text'; readonly text: string }
  | { readonly kind: 'reference'; readonly reference: Reference };

export interface ToolDefinition {
  readonly tool: { readonly id: string };
  /** Human-facing label; the provider-safe alias is derived from `tool.id`. */
  readonly name: string;
  readonly description: string;
  readonly inputSchema: JsonObject;
}

export interface ToolCallPart {
  readonly kind: 'tool-call';
  readonly requestId: string;
  readonly tool: { readonly id: string };
  readonly arguments: Json;
  /** Provider replay data (signatures, encrypted reasoning), kept with the message. */
  readonly providerOptions?: JsonObject;
}

export interface ToolResultPart {
  readonly kind: 'tool-result';
  readonly requestId: string;
  readonly tool: { readonly id: string };
  readonly status: 'succeeded' | 'failed';
  readonly parts: readonly ModelVisiblePart[];
}

/** Declares a tool; with `deferred: false` it also loads it. Adapters never send it as text. */
export interface ToolDefinitionPart {
  readonly kind: 'tool-definition';
  readonly definition: ToolDefinition;
  readonly deferred: boolean;
}

export interface SystemMessage {
  readonly role: 'system';
  readonly parts: readonly (ModelVisiblePart | ToolDefinitionPart)[];
}

export interface UserMessage {
  readonly role: 'user';
  readonly parts: readonly (ModelVisiblePart | ToolResultPart)[];
}

export interface AssistantMessage {
  readonly role: 'assistant';
  readonly parts: readonly (ModelVisiblePart | ToolCallPart)[];
}

export type Message = SystemMessage | UserMessage | AssistantMessage;
