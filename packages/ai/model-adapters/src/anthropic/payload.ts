/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { JsonObject, ToolDefinition } from '@molstar/ai-agent';
import { UnrepresentableInput } from '@molstar/ai-agent';
import {
  type ContextItem,
  encodeMessage,
  encodeVisible,
  estimateTool,
  type MeasurementPolicy,
  type NativeProjection,
  type ProjectionPolicy,
  uniqueTools,
  visibleParts,
  wireName,
} from '../common/payload.js';
import { projectVisibleContent } from './attachments.js';
import { assertSupportedParameters, supportsAnthropicAppendOnly } from './capabilities.js';

/** SDK-native reference loading. The inventory is the tool list, deferred except what loads before the
 * conversation starts; a later load becomes a `tool_addition` at its position. */
async function nativeAnthropicContext(
  items: readonly ContextItem[],
  inventory: readonly ToolDefinition[],
  model: string,
): Promise<NativeProjection> {
  const initial = new Set<string>();
  const system: string[] = [];
  const messages: JsonObject[] = [];
  let started = false;
  let inlineText: string[] = [];
  let changes: JsonObject[] = [];

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
      text: inlineText.join('\n\n'),
      ...(changes.length === 0 ? {} : { toolChanges: changes }),
    });
    inlineText = [];
    changes = [];
  }

  for (const item of items) {
    if (item.kind === 'load') {
      const id = item.definition.tool.id;
      if (!started) initial.add(id);
      else changes.push({ type: 'tool_addition', toolName: wireName(id) });
      continue;
    }
    if (item.message.role === 'system') {
      const encoded = await encodeVisible(visibleParts(item.message), model, 'system', projectVisibleContent);
      const text = encoded.text as string;
      if (text.length > 0) (started ? inlineText : system).push(text);
      continue;
    }
    // A tool load takes effect at the end of its user turn (later user messages in that turn go
    // first); in-place system text stays where it is written.
    if (item.message.role === 'assistant' || inlineText.length > 0) flush(item.message.role);
    started = true;
    messages.push(...(await encodeMessage(item.message, model, projectVisibleContent)));
  }
  flush();
  return {
    system: system.join('\n\n'),
    messages,
    tools: uniqueTools(inventory).map((tool) => ({ ...tool, ...(initial.has(tool.id) ? {} : { deferLoading: true }) })),
  };
}

export const projectionPolicy = (cache: boolean): ProjectionPolicy => ({
  visible: projectVisibleContent,
  assertParameters: assertSupportedParameters,
  native: { supports: supportsAnthropicAppendOnly, project: nativeAnthropicContext },
  fields: { cache },
});

export const measurementPolicy: MeasurementPolicy = {
  // Claude injects tool-use instructions: its published one-tool count_tokens example
  // records 403 tokens, beyond schema text alone. This fixed allowance is a heuristic.
  initialOverhead: (payload) => (payload.tools.length > 0 ? 313 : 0),
  disclosureTokens(message, tools) {
    return (Array.isArray(message.toolChanges) ? message.toolChanges : []).reduce<number>((sum, change) => {
      const name = (change as JsonObject).toolName;
      const tool = typeof name === 'string' ? tools.get(name) : undefined;
      if (tool === undefined) throw new Error('Projected disclosure references an undeclared tool.');
      return sum + estimateTool(tool);
    }, 0);
  },
};
