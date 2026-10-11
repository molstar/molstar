/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { JsonObject, ToolDefinition } from '@molstar/ai-agent';
import { UnrepresentableInput } from '@molstar/ai-agent';
import { estimateTextTokens } from '../common/estimate.js';
import {
  type ContextItem,
  encodeMessage,
  encodeVisible,
  type MeasurementPolicy,
  type NativeProjection,
  type ProjectionPolicy,
  uniqueTools,
  visibleParts,
} from '../common/payload.js';
import { projectVisibleContent } from './attachments.js';
import { assertSupportedParameters, supportsOpenAIAppendOnly } from './capabilities.js';

/** Full schemas render only at their load. Responses accepts positional full additions, so the tool
 * list the SDK gets serves parsing returned calls and is dropped before the request leaves. */
async function nativeOpenAIContext(
  items: readonly ContextItem[],
  inventory: readonly ToolDefinition[],
  model: string,
): Promise<NativeProjection> {
  const tools = uniqueTools(inventory);
  const pending = new Set<string>();
  const messages: JsonObject[] = [];
  for (const item of items) {
    if (item.kind === 'load') {
      if (pending.size > 0)
        throw new UnrepresentableInput('OpenAI tool disclosure requires a complete tool-result boundary.');
      const tool = uniqueTools([item.definition])[0];
      if (tool === undefined) continue;
      messages.push({
        type: 'additional_tools',
        role: 'developer',
        tools: [
          {
            type: 'function',
            name: tool.wireName,
            description: tool.description,
            parameters: tool.inputSchema,
            strict: false,
          },
        ],
      });
      continue;
    }
    const message = item.message;
    // Compaction may retain completed history while unloading its definitions: the active catalog,
    // not earlier calls, decides which tools can be called now.
    if (message.role === 'assistant') {
      for (const call of message.parts) if (call.kind === 'tool-call') pending.add(call.requestId);
    } else if (message.role === 'user') {
      for (const result of message.parts) if (result.kind === 'tool-result') pending.delete(result.requestId);
    } else if (pending.size > 0) {
      throw new UnrepresentableInput('OpenAI inline guidance requires a complete tool-result boundary.');
    }
    if (message.role === 'system') {
      messages.push({
        role: 'system',
        ...(await encodeVisible(visibleParts(message), model, 'system', projectVisibleContent)),
      });
    } else messages.push(...(await encodeMessage(message, model, projectVisibleContent)));
  }
  return { system: '', messages, tools };
}

export const projectionPolicy: ProjectionPolicy = {
  visible: projectVisibleContent,
  assertParameters: assertSupportedParameters,
  native: { supports: supportsOpenAIAppendOnly, project: nativeOpenAIContext },
  reasoning: (parameters) =>
    parameters.reasoning?.mode === 'off' ? { mode: 'adaptive', effort: 'none' } : parameters.reasoning,
};

export const measurementPolicy: MeasurementPolicy = {
  // The native catalog serves parsing only; prompt schemas are charged at their additions.
  includeTool: (_tool, payload) => payload.appendsTools !== true,
  disclosureTokens(message) {
    return (message.type === 'additional_tools' && Array.isArray(message.tools) ? message.tools : []).reduce<number>(
      (sum, addition) => {
        const tool = addition as JsonObject;
        return (
          sum +
          16 +
          estimateTextTokens(String(tool.name ?? '') + String(tool.description ?? '') + JSON.stringify(tool.parameters))
        );
      },
      0,
    );
  },
};
