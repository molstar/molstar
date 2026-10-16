/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  Json,
  JsonObject,
  Message,
  ModelVisiblePart,
  SystemMessage,
  ToolCallPart,
  ToolDefinition,
  ToolResultPart,
} from '../../messages/index.js';
import { canonicalJson, UnrepresentableInput } from '../../messages/index.js';
import { type AttachmentPlacement, type AttachmentRule, visibleParts } from './attachments.js';

// `Message` → the AI SDK's `ModelMessage`, once, at projection (providers.md §2). The result is JSON (the retained
// payload) and goes to `streamText` as it is. A message's label is never read.

/** Provider options as the SDK takes them: a record per provider. */
export type ProviderOptions = Readonly<Record<string, JsonObject>>;

/** A reasoning part, replayed before the call it led to. */
export interface ReasoningPart {
  readonly type: 'reasoning';
  readonly text: string;
  readonly providerOptions: ProviderOptions;
}

type Content = string | readonly JsonObject[];

/** A `ModelMessage` as JSON. */
export type SdkMessage =
  | { readonly role: 'system'; readonly content: string; readonly providerOptions?: ProviderOptions }
  | { readonly role: 'user' | 'assistant'; readonly content: Content }
  | { readonly role: 'tool'; readonly content: readonly JsonObject[] };

/** A tool as the request names it. */
export interface WireTool {
  readonly id: string;
  readonly name: string;
  readonly description: string;
  readonly inputSchema: JsonObject;
  readonly deferLoading?: true;
}

/** How a provider replays what a call carried (thinking signatures, encrypted reasoning, item ids). */
export interface ReplayCodec {
  /** The reasoning parts that go before a call, from its `providerOptions`. */
  reasoning?(options: Json | undefined): readonly ReasoningPart[];
  /** The tool-call part's own provider options. */
  callOptions?(options: Json | undefined): ProviderOptions | undefined;
}

export interface Converter {
  readonly model: string;
  readonly attachments: AttachmentRule;
  readonly replay: ReplayCodec;
}

export function wireName(toolId: string): string {
  const name = toolId.replace(/[^a-zA-Z0-9_-]/g, '_');
  if (name.length === 0) throw new Error('Tool id does not produce a provider-safe name.');
  return name;
}

export function uniqueTools(tools: readonly ToolDefinition[]): WireTool[] {
  const seen = new Set<string>();
  return tools.map((tool) => {
    const name = wireName(tool.tool.id);
    if (seen.has(name)) throw new Error(`Duplicate projected tool alias ${name}.`);
    seen.add(name);
    return { id: tool.tool.id, name, description: tool.description, inputSchema: tool.inputSchema };
  });
}

const visible = (c: Converter, parts: readonly ModelVisiblePart[], placement: AttachmentPlacement) =>
  visibleParts(parts, placement, c.model, c.attachments);

/** One text, or the parts when an attachment is among them. */
const contentOf = (v: ReturnType<typeof visible>): Content =>
  v.rich ? (v.parts as unknown as readonly JsonObject[]) : v.text;

/** A system message's text: its tool definitions are never prompt text. */
export function systemText(c: Converter, message: SystemMessage): string {
  const parts = message.parts.filter((part): part is ModelVisiblePart => part.kind !== 'tool-definition');
  return visible(c, parts, 'system').text;
}

function toolResult(c: Converter, part: ToolResultPart): JsonObject {
  const v = visible(c, part.parts, 'tool-result');
  const files = v.parts.some((p) => p.type === 'file');
  if (part.status === 'failed' && files) {
    throw new UnrepresentableInput('Binary attachments in failed tool-result content are unsupported.');
  }
  return {
    type: 'tool-result',
    toolCallId: part.requestId,
    toolName: wireName(part.tool.id),
    output:
      part.status === 'failed'
        ? { type: 'error-text', value: v.text }
        : v.rich
          ? { type: 'content', value: v.parts as unknown as Json[] }
          : { type: 'text', value: v.text },
  };
}

function toolCall(c: Converter, part: ToolCallPart): JsonObject {
  const options = c.replay.callOptions?.(part.providerOptions);
  return {
    type: 'tool-call',
    toolCallId: part.requestId,
    toolName: wireName(part.tool.id),
    input: part.arguments,
    ...(options === undefined ? {} : { providerOptions: options }),
  };
}

/** A user or assistant message as `ModelMessage`s (a system message is the projection's to place). Tool results
 * come first, in a tool message, and the user's other content after them: results must follow their calls. */
export function convertMessage(c: Converter, message: Message): SdkMessage[] {
  if (message.role === 'system') {
    return [{ role: 'system', content: systemText(c, message) }];
  }
  if (message.role === 'assistant') {
    const v = visible(
      c,
      message.parts.filter((part): part is ModelVisiblePart => part.kind !== 'tool-call'),
      'assistant',
    );
    const calls = message.parts.filter((part): part is ToolCallPart => part.kind === 'tool-call');
    if (calls.length === 0) return [{ role: 'assistant', content: contentOf(v) }];
    const reasoning = (call: ToolCallPart) =>
      (c.replay.reasoning?.(call.providerOptions) ?? []) as unknown as JsonObject[];
    const text = v.rich
      ? (v.parts as unknown as JsonObject[])
      : v.text.length > 0
        ? [{ type: 'text', text: v.text }]
        : [];
    return [
      {
        role: 'assistant',
        content: [
          // Provider reasoning precedes the call it accompanied.
          ...reasoning(calls[0] as ToolCallPart),
          ...text,
          ...calls.flatMap((call, i) => [...(i > 0 ? reasoning(call) : []), toolCall(c, call)]),
        ],
      },
    ];
  }
  const results = message.parts.filter((part): part is ToolResultPart => part.kind === 'tool-result');
  const v = visible(
    c,
    message.parts.filter((part): part is ModelVisiblePart => part.kind !== 'tool-result'),
    'user',
  );
  if (results.length === 0) return [{ role: 'user', content: contentOf(v) }];
  return [
    { role: 'tool', content: results.map((result) => toolResult(c, result)) },
    ...(v.text !== '' || v.rich ? [{ role: 'user' as const, content: contentOf(v) }] : []),
  ];
}

/** The collected layout: system text only before the conversation, every available tool in the list. */
export function collected(
  c: Converter & { readonly provider: string },
  messages: readonly Message[],
  tools: readonly ToolDefinition[],
): { system: string; messages: SdkMessage[]; tools: WireTool[] } {
  const system: string[] = [];
  const out: SdkMessage[] = [];
  for (const message of messages) {
    if (message.role !== 'system') {
      out.push(...convertMessage(c, message));
      continue;
    }
    const text = systemText(c, message);
    if (text.length === 0) continue;
    if (out.length > 0) {
      throw new UnrepresentableInput(
        `${c.provider} model ${c.model} cannot place system text inside the conversation.`,
      );
    }
    system.push(text);
  }
  return { system: system.join('\n\n'), messages: out, tools: uniqueTools(tools) };
}

/** The system prompt, loads and messages in the order a head holds them. */
export type ContextItem =
  | { readonly kind: 'message'; readonly message: Message }
  | { readonly kind: 'load'; readonly definition: ToolDefinition };

const same = (a: ToolDefinition, b: ToolDefinition): boolean =>
  canonicalJson(a as unknown as Json) === canonicalJson(b as unknown as Json);

/** The messages as ordered context, for the appended layout: a loading `tool-definition` part becomes a `load` at its
 * position (tools given without a part load first), and deferred parts only declare. Loads follow the head, not what
 * this request may call, so a request with fewer tools (a nested run) keeps the cached prefix; a call to a tool it may
 * not call is unbound (harness.md §8). `inventory` is every definition the head names, in order of appearance. */
export function orderContext(
  messages: readonly Message[],
  tools: readonly ToolDefinition[],
): { readonly items: readonly ContextItem[]; readonly inventory: readonly ToolDefinition[] } {
  const inventory = new Map<string, ToolDefinition>();
  const loaded = new Set<string>();
  const items: ContextItem[] = [];
  for (const message of messages) {
    if (message.role !== 'system') {
      items.push({ kind: 'message', message });
      continue;
    }
    let run: ModelVisiblePart[] = [];
    const flush = () => {
      if (run.length > 0) items.push({ kind: 'message', message: { role: 'system', parts: run } });
      run = [];
    };
    for (const part of message.parts) {
      if (part.kind !== 'tool-definition') {
        run.push(part);
        continue;
      }
      const { definition, deferred } = part;
      const known = inventory.get(definition.tool.id);
      if (known === undefined) inventory.set(definition.tool.id, definition);
      else if (!same(known, definition)) {
        throw new Error(
          `Native tool references cannot represent multiple definition versions of ${definition.tool.id}; use collected projection or a new context head.`,
        );
      }
      if (deferred || loaded.has(definition.tool.id)) continue;
      flush();
      loaded.add(definition.tool.id);
      items.push({ kind: 'load', definition });
    }
    flush();
  }
  const initial = tools.filter((tool) => !loaded.has(tool.tool.id));
  for (const tool of initial) inventory.set(tool.tool.id, tool);
  return {
    items: [...initial.map((definition): ContextItem => ({ kind: 'load', definition })), ...items],
    inventory: [...inventory.values()],
  };
}
