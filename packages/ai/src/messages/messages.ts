/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, JsonObject, Reference } from './values.js';

/** What a model sees of content: text, or a reference the provider resolves. */
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

/** Declares a tool; with `deferred: false` it also loads it. Providers never send it as text. */
export interface ToolDefinitionPart {
  readonly kind: 'tool-definition';
  readonly definition: ToolDefinition;
  readonly deferred: boolean;
}

/**
 * What an entry is, for the code that reads the journal (a session, compaction, a trace): its request, its guidance, a
 * skill's guidance. Never sent to a provider. Part of the message, so of its hash: it says only what is true of the
 * content, never of an occurrence (a turn, a time, a step, a cost: those are the event log's, harness.md §1).
 */
export interface EntryLabel {
  /** Namespaced by the package that writes it: `session:request`, `tools:skill`. */
  readonly kind: string;
  /** Which one within the kind: a skill's `skill:animation@3f2a9c`. */
  readonly key?: string;
  /** Small facts the parts also say, as data: a request's prompt and focus. */
  readonly data?: JsonObject;
}

export interface SystemMessage {
  readonly role: 'system';
  readonly parts: readonly (ModelVisiblePart | ToolDefinitionPart)[];
  readonly label?: EntryLabel;
}

export interface UserMessage {
  readonly role: 'user';
  readonly parts: readonly (ModelVisiblePart | ToolResultPart)[];
  readonly label?: EntryLabel;
}

export interface AssistantMessage {
  readonly role: 'assistant';
  readonly parts: readonly (ModelVisiblePart | ToolCallPart)[];
  readonly label?: EntryLabel;
}

export type Message = SystemMessage | UserMessage | AssistantMessage;
