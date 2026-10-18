/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, JsonObject, ToolDefinition } from '../../messages/index.js';
import { UnrepresentableInput } from '../../messages/index.js';
import { estimateTool } from '../common/estimate.js';
import {
  type ContextItem,
  type Converter,
  convertMessage,
  type SdkMessage,
  systemText,
  uniqueTools,
} from '../common/messages.js';
import type { Layout, Payload } from '../common/provider.js';

// OpenAI's positional tool additions (`additional_tools` developer items), for the catalog's models with the
// `appendsTools` quirk (`_spec/models.md`, "Projection rules"). The installed SDK has no `additional_tools` input: a load is a
// placeholder system message (its tools in `providerOptions.molstar`), which `additionsTransport` replaces in the
// request body before it leaves.

/** The placeholder for the `n`-th addition of a request. */
const marker = (n: number) => `molstar-ai:openai:additional-tools:${n}`;

type Addition = { readonly type: 'function'; readonly name: string } & JsonObject;

/** The tool list is the head's inventory and serves only to parse returned calls (the transport drops it); a full
 * schema renders where it loads. */
export function appended(c: Converter, items: readonly ContextItem[], inventory: readonly ToolDefinition[]): Layout {
  const pending = new Set<string>();
  const messages: SdkMessage[] = [];
  let additions = 0;
  for (const item of items) {
    if (item.kind === 'load') {
      if (pending.size > 0)
        throw new UnrepresentableInput('OpenAI tool disclosure requires a complete tool-result boundary.');
      const [tool] = uniqueTools([item.definition]);
      if (tool === undefined) continue;
      const addition: Addition = {
        type: 'function',
        name: tool.name,
        description: tool.description,
        parameters: tool.inputSchema,
        strict: false,
      };
      messages.push({
        role: 'system',
        content: marker(additions++),
        providerOptions: { molstar: { additionalTools: [addition] } },
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
    if (message.role === 'system') messages.push({ role: 'system', content: systemText(c, message) });
    else messages.push(...convertMessage(c, message));
  }
  return { system: '', messages, tools: uniqueTools(inventory) };
}

const additionsOf = (message: SdkMessage): readonly Addition[] => {
  if (message.role !== 'system') return [];
  const tools = (message.providerOptions?.molstar as { additionalTools?: Json } | undefined)?.additionalTools;
  return Array.isArray(tools) ? (tools as unknown as Addition[]) : [];
};

/** An addition's schema, where it loads. */
export const additionTokens = (message: SdkMessage): number =>
  additionsOf(message).reduce(
    (sum, tool) =>
      sum +
      estimateTool({
        name: String(tool.name ?? ''),
        description: String(tool.description ?? ''),
        inputSchema: tool.parameters,
      }),
    0,
  );

/** The placeholders become `additional_tools` items, in order, and the tool list leaves the body. */
export function additionsTransport(payload: Payload, fetch: typeof globalThis.fetch): typeof globalThis.fetch {
  const additions = payload.messages.flatMap((message) => {
    const tools = additionsOf(message);
    return tools.length > 0 ? [{ marker: message.content as string, tools }] : [];
  });
  if (payload.appended !== true) return fetch;
  return async (url, init) => {
    if (typeof init?.body !== 'string') throw new Error('OpenAI disclosure transport expected an SDK JSON request.');
    const body = JSON.parse(init.body) as Record<string, unknown>;
    if (!Array.isArray(body.input)) throw new Error('OpenAI disclosure transport expected ordered input items.');
    let next = 0;
    body.input = body.input.map((entry) => {
      if (entry === null || typeof entry !== 'object' || Array.isArray(entry)) return entry;
      const item = entry as JsonObject;
      const candidate = additions.find((addition) => addition.marker === item.content);
      if (candidate === undefined) return entry;
      if (item.role !== 'developer' || candidate !== additions[next])
        throw new Error('OpenAI disclosure placeholder order or role changed during SDK conversion.');
      next++;
      return { type: 'additional_tools', role: 'developer', tools: candidate.tools };
    });
    if (next !== additions.length) throw new Error('OpenAI disclosure placeholder was lost during SDK conversion.');
    // The tools were given to the SDK only to parse returned calls; their full definitions stand at their positions.
    delete body.tools;
    return fetch(url, { ...init, body: JSON.stringify(body) });
  };
}
