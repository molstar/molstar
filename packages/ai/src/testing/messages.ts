/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AssistantMessage, Json, Message, SystemMessage, ToolDefinition, UserMessage } from '../messages/index.js';

export const user = (text: string): UserMessage => ({ role: 'user', parts: [{ kind: 'text', text }] });
export const system = (text: string): SystemMessage => ({ role: 'system', parts: [{ kind: 'text', text }] });
export const assistant = (text: string): AssistantMessage => ({ role: 'assistant', parts: [{ kind: 'text', text }] });

/** A tool definition with an object input schema. */
export function toolDefinition(id: string, properties: Record<string, Json> = {}): ToolDefinition {
  return {
    tool: { id },
    name: id,
    description: `The ${id} tool.`,
    inputSchema: { type: 'object', properties },
  };
}

/** A system message that declares (and, unless deferred, loads) the given tools. */
export function declareTools(definitions: readonly ToolDefinition[], deferred = false): SystemMessage {
  return {
    role: 'system',
    parts: definitions.map((definition) => ({ kind: 'tool-definition', definition, deferred })),
  };
}

/** The text of a message's text parts, joined. */
export function textOf(message: Message): string {
  return message.parts
    .map((part) => (part.kind === 'text' ? part.text : ''))
    .filter((text) => text !== '')
    .join('\n');
}
