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
  ModelParameters,
  ModelVisiblePart,
  ProjectInput,
  ReasoningParameters,
  SystemMessage,
  ToolCallPart,
  ToolDefinition,
  ToolResultPart,
} from '../../messages/index.js';
import { canonicalJson, UnrepresentableInput } from '../../messages/index.js';
import { estimateAttachmentTokens, type VisibleProjector } from './attachments.js';
import { defaultMaxOutputTokens, validateCapabilities } from './capabilities.js';
import { estimateTextTokens, type TokenEstimate } from './estimate.js';

export interface WireTool {
  readonly id: string;
  readonly wireName: string;
  readonly description: string;
  readonly inputSchema: JsonObject;
  readonly deferLoading?: boolean;
}

/** What an adapter retains of a request: neutral, JSON, and never the API key. */
export interface ProjectedPayload {
  readonly adapterRevision: string;
  readonly provider: string;
  readonly model: string;
  readonly system: string;
  readonly messages: readonly JsonObject[];
  readonly tools: readonly WireTool[];
  readonly maxOutputTokens: number;
  readonly reasoning?: Exclude<ReasoningParameters, { readonly mode: 'off' }>;
  readonly temperature?: number;
  readonly topP?: number;
  readonly cache?: boolean;
  /** Only where the capabilities map it (OpenAI's prompt cache key). */
  readonly cacheScopeId?: string;
  /** Tool definitions load where their parts stand, in the provider's own way. */
  readonly appendsTools?: true;
}

export function wireName(toolId: string): string {
  const name = toolId.replace(/[^a-zA-Z0-9_-]/g, '_');
  if (name.length === 0) throw new Error('Tool id does not produce a provider-safe name.');
  return name;
}

/** The system prompt, loads and messages in the order a head holds them. */
export type ContextItem =
  | { readonly kind: 'message'; readonly message: Message }
  | { readonly kind: 'load'; readonly definition: ToolDefinition };

export interface NativeProjection {
  readonly system: string;
  readonly messages: readonly JsonObject[];
  readonly tools: readonly WireTool[];
}

export interface ProjectionPolicy {
  readonly visible: VisibleProjector;
  readonly assertParameters: (parameters: ModelParameters, model: string) => void;
  readonly reasoning?: (parameters: ModelParameters) => ReasoningParameters | undefined;
  /** The provider's way of adding a tool mid-conversation, for the models that have it. */
  readonly native?: {
    readonly supports: (model: string) => boolean;
    readonly project: (
      items: readonly ContextItem[],
      inventory: readonly ToolDefinition[],
      model: string,
    ) => Promise<NativeProjection>;
  };
  readonly fields?: Pick<ProjectedPayload, 'cache'>;
}

const same = (a: ToolDefinition, b: ToolDefinition): boolean =>
  canonicalJson(a as unknown as Json) === canonicalJson(b as unknown as Json);

/** The messages as ordered context: a loading `tool-definition` part becomes a `load` at its position
 * (tools given without a part load first), and deferred parts only declare. Loads follow the head,
 * not what this request may call, so a request with fewer tools (a nested run) keeps the cached
 * prefix; a call to a tool it may not call is unbound (harness.md §8). `inventory` is every definition
 * the head names, in order of appearance. */
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

export async function buildProjectedPayload(
  input: ProjectInput,
  provider: string,
  adapterRevision: string,
  policy: ProjectionPolicy,
): Promise<ProjectedPayload> {
  const { configuration, capabilities } = input;
  const model = configuration.target.model;
  if (configuration.target.provider !== provider) {
    throw new Error(`${provider} adapter cannot project provider ${configuration.target.provider}.`);
  }
  validateCapabilities(capabilities, configuration.parameters);
  policy.assertParameters(configuration.parameters, model);
  const reasoning = policy.reasoning?.(configuration.parameters) ?? configuration.parameters.reasoning;
  let layout: NativeProjection;
  if (capabilities.appendsTools) {
    if (policy.native?.supports(model) !== true) {
      throw new Error(`${provider} model ${model} does not support appended tool definitions.`);
    }
    const { items, inventory } = orderContext(input.messages, input.tools);
    layout = await policy.native.project(items, inventory, model);
  } else {
    const system: string[] = [];
    const messages: JsonObject[] = [];
    for (const message of input.messages) {
      if (message.role !== 'system') {
        messages.push(...(await encodeMessage(message, model, policy.visible)));
        continue;
      }
      const text = contentText(await policy.visible(visibleParts(message), model, 'system'));
      if (text.length === 0) continue;
      if (messages.length > 0) {
        throw new UnrepresentableInput(`${provider} model ${model} cannot place system text inside the conversation.`);
      }
      system.push(text);
    }
    layout = { system: system.join('\n\n'), messages, tools: uniqueTools(input.tools) };
  }
  const { parameters } = configuration;
  return {
    adapterRevision,
    provider,
    model,
    ...layout,
    maxOutputTokens: parameters.maxOutputTokens ?? capabilities.defaultMaxOutputTokens ?? defaultMaxOutputTokens,
    ...(reasoning !== undefined && reasoning.mode !== 'off' ? { reasoning } : {}),
    ...(parameters.temperature !== undefined ? { temperature: parameters.temperature } : {}),
    ...(parameters.topP !== undefined ? { topP: parameters.topP } : {}),
    ...policy.fields,
    ...(capabilities.cacheScope && configuration.cacheScopeId !== undefined
      ? { cacheScopeId: configuration.cacheScopeId }
      : {}),
    ...(capabilities.appendsTools ? { appendsTools: true as const } : {}),
  };
}

/** A retained payload, checked for the fields invocation reads. */
export function decodeProjectedPayload(payload: Json): ProjectedPayload {
  const value = payload as { [key: string]: Json } | null;
  if (
    value === null ||
    typeof value !== 'object' ||
    Array.isArray(value) ||
    typeof value.provider !== 'string' ||
    typeof value.model !== 'string' ||
    typeof value.adapterRevision !== 'string' ||
    typeof value.maxOutputTokens !== 'number' ||
    !Array.isArray(value.messages) ||
    !Array.isArray(value.tools)
  ) {
    throw new Error('Projected payload is malformed.');
  }
  return { system: '', ...value } as unknown as ProjectedPayload;
}

export interface MeasurementPolicy {
  readonly initialOverhead?: (payload: ProjectedPayload) => number;
  readonly disclosureTokens?: (message: JsonObject, tools: ReadonlyMap<string, WireTool>) => number;
  readonly includeTool?: (tool: WireTool, payload: ProjectedPayload) => boolean;
  /** The application's calibration for this provider (default: four characters a token). */
  readonly profile?: TokenEstimate;
}

/** The input estimate for admission: active prompt content, not transport JSON or output controls.
 * Per-message and tool overheads are approximate and deterministic. */
export function measureProjected(payload: ProjectedPayload, policy: MeasurementPolicy = {}): number {
  const toolsByName = new Map(payload.tools.map((tool) => [tool.wireName, tool]));
  const text = (value: string) => estimateTextTokens(value, policy.profile);
  let tokens = text(payload.system) + (policy.initialOverhead?.(payload) ?? 0);
  for (const message of payload.messages) {
    tokens += 8 + estimateVisible(message, text);
    for (const call of Array.isArray(message.toolCalls) ? message.toolCalls : []) {
      const item = call as JsonObject;
      tokens += 8 + text(JSON.stringify(item.arguments ?? {}));
      if (typeof item.wireName === 'string') tokens += text(item.wireName);
      if (item.providerOptions !== undefined) tokens += estimateThinkingText(item.providerOptions, text);
    }
    for (const result of Array.isArray(message.toolResults) ? message.toolResults : []) {
      tokens += 8 + estimateVisible(result as JsonObject, text);
    }
    tokens += policy.disclosureTokens?.(message, toolsByName) ?? 0;
  }
  for (const tool of payload.tools) {
    if (policy.includeTool?.(tool, payload) ?? tool.deferLoading !== true) tokens += estimateTool(tool, policy.profile);
  }
  return Math.max(1, tokens);
}

export function estimateTool(tool: WireTool, profile?: TokenEstimate): number {
  return 16 + estimateTextTokens(tool.wireName + tool.description + JSON.stringify(tool.inputSchema), profile);
}

function estimateVisible(value: JsonObject, text: (value: string) => number): number {
  if (!Array.isArray(value.content)) return text(typeof value.text === 'string' ? value.text : '');
  return value.content.reduce<number>((tokens, item) => {
    const part = item as JsonObject;
    return tokens + 4 + text(typeof part.text === 'string' ? part.text : '') + estimateAttachmentTokens(part);
  }, 0);
}

/** Replayed thinking text counts; opaque signatures do not. */
function estimateThinkingText(value: Json, text: (value: string) => number): number {
  if (Array.isArray(value)) return value.reduce<number>((total, item) => total + estimateThinkingText(item, text), 0);
  if (value === null || typeof value !== 'object') return 0;
  return Object.entries(value).reduce(
    (total, [key, item]) =>
      total +
      ((key === 'thinking' || key === 'text') && typeof item === 'string'
        ? text(item)
        : estimateThinkingText(item, text)),
    0,
  );
}

export function uniqueTools(tools: readonly ToolDefinition[]): readonly WireTool[] {
  const seen = new Set<string>();
  return tools.map((tool) => {
    const name = wireName(tool.tool.id);
    if (seen.has(name)) throw new Error(`Duplicate projected tool alias ${name}.`);
    seen.add(name);
    return { id: tool.tool.id, wireName: name, description: tool.description, inputSchema: tool.inputSchema };
  });
}

/** A system message's text and attachments: its tool definitions are never prompt text. */
export function visibleParts(message: SystemMessage): readonly ModelVisiblePart[] {
  return message.parts.filter((part): part is ModelVisiblePart => part.kind !== 'tool-definition');
}

export async function encodeMessage(
  message: Message,
  model: string,
  visible: VisibleProjector,
): Promise<readonly JsonObject[]> {
  if (message.role === 'assistant') {
    return [
      {
        role: 'assistant',
        ...(await encodeVisible(
          message.parts.filter((part): part is ModelVisiblePart => part.kind !== 'tool-call'),
          model,
          'assistant',
          visible,
        )),
        toolCalls: message.parts.filter((part): part is ToolCallPart => part.kind === 'tool-call').map(encodeToolCall),
      },
    ];
  }
  if (message.role === 'user') {
    const results = message.parts.filter((part): part is ToolResultPart => part.kind === 'tool-result');
    const content = await encodeVisible(
      message.parts.filter((part): part is ModelVisiblePart => part.kind !== 'tool-result'),
      model,
      'user',
      visible,
    );
    if (results.length === 0) return [{ role: 'user', ...content }];
    // Preserve required result adjacency before any visible follow-up content.
    const encoded = [];
    for (const result of results) encoded.push(await encodeToolResult(result, model, visible));
    return [
      { role: 'tool', toolResults: encoded },
      ...(content.text !== '' || content.content !== undefined ? [{ role: 'user', ...content }] : []),
    ];
  }
  return [{ role: 'system', ...(await encodeVisible(visibleParts(message), model, 'system', visible)) }];
}

export function contentText(content: readonly JsonObject[]): string {
  return content.map((part) => (typeof part.text === 'string' ? part.text : '')).join('');
}

export async function encodeVisible(
  parts: readonly ModelVisiblePart[],
  model: string,
  placement: 'system' | 'user' | 'assistant' | 'tool-result',
  visible: VisibleProjector,
): Promise<JsonObject> {
  const content = await visible(parts, model, placement);
  // Keep the text-only representation, and retain ordered blocks for files.
  return content.some((part) => part.type === 'file' || part.attachment !== undefined)
    ? { content, text: contentText(content) }
    : { text: contentText(content) };
}

function encodeToolCall(part: ToolCallPart): JsonObject {
  return {
    requestId: part.requestId,
    toolId: part.tool.id,
    wireName: wireName(part.tool.id),
    arguments: part.arguments,
    ...(part.providerOptions === undefined ? {} : { providerOptions: part.providerOptions }),
  };
}

async function encodeToolResult(part: ToolResultPart, model: string, visible: VisibleProjector): Promise<JsonObject> {
  const content = await encodeVisible(part.parts, model, 'tool-result', visible);
  if (
    part.status === 'failed' &&
    Array.isArray(content.content) &&
    content.content.some((item) => (item as JsonObject).type === 'file')
  ) {
    throw new UnrepresentableInput('Binary attachments in failed tool-result content are unsupported.');
  }
  return {
    requestId: part.requestId,
    toolId: part.tool.id,
    wireName: wireName(part.tool.id),
    status: part.status,
    ...content,
  };
}
