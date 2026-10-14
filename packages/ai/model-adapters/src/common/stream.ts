/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  AssistantMessage,
  Json,
  JsonObject,
  ModelFinish,
  ModelFrame,
  ProjectedRequest,
  ProviderUsage,
  ToolCallPart,
} from '@molstar/ai-agent';
import { digest, ModelError } from '@molstar/ai-agent';
import { jsonSchema, type LanguageModelUsage, type ModelMessage, streamText, type TextPart, type ToolSet } from 'ai';
import { decodeProjectedPayload, type ProjectedPayload, type WireTool } from './payload.js';

export interface ReplayState {
  reasoning(part: StreamPart): void;
  toolCall(call: ToolCallPart, part: StreamPart, index: number): ToolCallPart;
  finish(calls: readonly ToolCallPart[]): readonly ToolCallPart[];
}
export interface StreamPolicy {
  options(payload: ProjectedPayload): Pick<Parameters<typeof streamText>[0], 'system' | 'providerOptions'>;
  toolOptions?(tool: WireTool, index: number, payload: ProjectedPayload): Record<string, JsonObject> | undefined;
  inlineMessage(message: JsonObject, index: number, digest: string): ModelMessage | undefined;
  reasoning?(value: Json | undefined): ReasoningPart[];
  replayOptions?(value: Json | undefined): Record<string, JsonObject> | undefined;
  createReplay(): ReplayState;
}

/** A failure before any request: the request is refused or malformed, and nothing was generated. */
const refused = (message: string) => new ModelError(message, { transient: false, effect: 'none' });

export async function* invokeProjectedRequest(input: {
  readonly request: ProjectedRequest;
  readonly signal: AbortSignal;
  readonly expectedRevision: string;
  readonly expectedProvider: string;
  readonly policy: StreamPolicy;
  readonly languageModel: Parameters<typeof streamText>[0]['model'];
}): AsyncIterable<ModelFrame> {
  if (input.request.adapterRevision !== input.expectedRevision) {
    throw refused('Projected request belongs to an incompatible adapter revision.');
  }
  if (input.request.digest !== (await digest(input.request.payload))) {
    throw refused('Projected request payload digest does not match the payload.');
  }
  const payload = decodeProjectedPayload(input.request.payload);
  if (payload.provider !== input.expectedProvider || payload.adapterRevision !== input.expectedRevision) {
    throw refused('Projected request provider does not match this adapter.');
  }
  const tools = sdkTools(payload, input.policy);
  const messages = sdkMessages(payload, input.request.digest, input.policy);
  let streamError: unknown;
  const streamOptions: Parameters<typeof streamText>[0] = {
    model: input.languageModel,
    abortSignal: input.signal,
    maxRetries: 0,
    maxOutputTokens: payload.maxOutputTokens,
    messages,
    ...(messages.some((message) => message.role === 'system') ? { allowSystemInMessages: true } : {}),
    onError({ error }) {
      streamError = error;
    },
    ...input.policy.options(payload),
    ...(payload.tools.length > 0 ? { tools } : {}),
    ...(payload.temperature !== undefined ? { temperature: payload.temperature } : {}),
    ...(payload.topP !== undefined ? { topP: payload.topP } : {}),
  };
  const result = streamText(streamOptions);
  const text: string[] = [];
  const toolCalls: ToolCallPart[] = [];
  let finish: ModelFinish | undefined;
  let usage: ProviderUsage = {};
  const replay = input.policy.createReplay();
  const parts = result.fullStream[Symbol.asyncIterator]();
  const aborted = waitForAbort(input.signal);
  try {
    while (true) {
      const next = await Promise.race([parts.next(), aborted]);
      if (next.done === true) break;
      const typed = next.value as StreamPart;
      if (typed.type === 'error') throw streamFailure(typed.error);
      if (typed.type === 'abort') throw streamFailure(typed.reason ?? 'aborted');
      if (typed.type === 'text-delta') {
        const delta = streamTextDelta(typed);
        if (delta.length > 0) {
          text.push(delta);
          yield { kind: 'text-delta', text: delta };
        }
        continue;
      }
      if (['reasoning-start', 'reasoning-delta', 'reasoning-end'].includes(typed.type)) replay.reasoning(typed);
      if (typed.type === 'reasoning-delta') {
        const delta = streamTextDelta(typed);
        if (delta.length > 0) yield { kind: 'reasoning-delta', text: delta };
        continue;
      }
      if (typed.type === 'tool-call') {
        toolCalls.push(replay.toolCall(streamToolCall(typed, payload.tools), typed, toolCalls.length));
        continue;
      }
      if (typed.type === 'finish') {
        finish = mapFinishReason(typed.finishReason);
        const nextUsage = streamUsage(typed);
        if (nextUsage !== undefined) {
          usage = nextUsage;
          yield { kind: 'usage', usage };
        }
      }
    }
  } finally {
    await parts.return?.();
  }
  if (streamError !== undefined) throw streamFailure(streamError);
  if (finish === undefined) throw new Error('Model stream ended without a recognized finish.');
  let finalUsage: LanguageModelUsage | undefined;
  try {
    finalUsage = await result.usage;
  } catch {
    finalUsage = undefined;
  }
  if (finalUsage !== undefined) usage = { ...usage, ...providerUsageFromSdk(finalUsage) };
  const message: AssistantMessage = {
    role: 'assistant',
    parts: [
      ...(text.length > 0 ? ([{ kind: 'text' as const, text: text.join('') }] as const) : []),
      ...replay.finish(toolCalls),
    ],
  };
  yield { kind: 'result', message, finish, usage };
}

function sdkTools(payload: ProjectedPayload, policy: StreamPolicy): ToolSet {
  const projected: ToolSet = {};
  for (const [index, tool] of payload.tools.entries()) {
    const providerOptions = policy.toolOptions?.(tool, index, payload);
    projected[tool.wireName] = {
      description: tool.description,
      inputSchema: jsonSchema(tool.inputSchema),
      ...(providerOptions === undefined ? {} : { providerOptions }),
    };
  }
  return projected;
}

function sdkMessages(payload: ProjectedPayload, requestDigest: string, policy: StreamPolicy): ModelMessage[] {
  const messages: ModelMessage[] = [];
  for (const [messageIndex, message] of payload.messages.entries()) {
    const role = message.role;
    const text = typeof message.text === 'string' ? message.text : '';
    const visible = sdkVisibleContent(message);
    const inline = policy.inlineMessage(message, messageIndex, requestDigest);
    if (inline !== undefined) {
      messages.push(inline);
      continue;
    }
    if (role === 'assistant') {
      const toolCalls = Array.isArray(message.toolCalls) ? message.toolCalls : [];
      if (toolCalls.length === 0) {
        messages.push({ role: 'assistant', content: Array.isArray(message.content) ? visible : text });
        continue;
      }
      // Provider reasoning must precede the function call it accompanied.
      const reasoning = (call: Json): ReasoningPart[] => policy.reasoning?.((call as JsonObject).providerOptions) ?? [];
      messages.push({
        role: 'assistant',
        content: [
          ...reasoning(toolCalls[0] ?? {}),
          ...visible,
          ...toolCalls.flatMap((call, i) => {
            const record = call as JsonObject;
            const providerOptions = policy.replayOptions?.(record.providerOptions);
            return [
              ...(i > 0 ? reasoning(call) : []),
              {
                type: 'tool-call' as const,
                toolCallId: String(record.requestId ?? ''),
                toolName: String(record.wireName ?? ''),
                input: record.arguments,
                ...(providerOptions === undefined ? {} : { providerOptions }),
              },
            ];
          }),
        ],
      });
      continue;
    }
    if (role === 'tool') {
      const results = Array.isArray(message.toolResults) ? message.toolResults : [];
      messages.push({
        role: 'tool',
        content: results.map((item) => {
          const record = item as JsonObject;
          const value = typeof record.text === 'string' ? record.text : '';
          return {
            type: 'tool-result' as const,
            toolCallId: String(record.requestId ?? ''),
            toolName: String(record.wireName ?? ''),
            output:
              record.status === 'failed'
                ? { type: 'error-text' as const, value }
                : Array.isArray(record.content)
                  ? { type: 'content' as const, value: sdkVisibleContent(record) }
                  : { type: 'text' as const, value },
          };
        }),
      });
      if (text.length > 0) messages.push({ role: 'user', content: text });
      continue;
    }
    messages.push({ role: 'user', content: Array.isArray(message.content) ? visible : text });
  }
  return messages;
}

type SdkVisiblePart =
  | TextPart
  | {
      type: 'file';
      data: { type: 'data'; data: string } | { type: 'text'; text: string };
      mediaType: string;
      filename?: string;
    };

function sdkVisibleContent(message: JsonObject): SdkVisiblePart[] {
  if (!Array.isArray(message.content)) {
    return typeof message.text === 'string' && message.text.length > 0 ? [{ type: 'text', text: message.text }] : [];
  }
  return message.content.map((item): SdkVisiblePart => {
    const part = item as JsonObject;
    if (part.type === 'text' && typeof part.text === 'string') return { type: 'text', text: part.text };
    if (part.type !== 'file' || typeof part.data !== 'string' || typeof part.mediaType !== 'string') {
      throw new Error('Invalid projected attachment part.');
    }
    return {
      type: 'file',
      data: typeof part.text === 'string' ? { type: 'text', text: part.text } : { type: 'data', data: part.data },
      mediaType: part.mediaType,
      ...(typeof part.name === 'string' ? { filename: part.name } : {}),
    };
  });
}

export interface StreamPart {
  readonly type: string;
  readonly text?: string;
  readonly delta?: string;
  readonly toolCallId?: string;
  readonly toolName?: string;
  readonly input?: unknown;
  readonly finishReason?: string;
  readonly error?: unknown;
  readonly reason?: unknown;
  readonly totalUsage?: LanguageModelUsage;
  readonly id?: string;
  readonly providerMetadata?: Readonly<Record<string, Readonly<Record<string, unknown>>>>;
}

/** An assistant message's reasoning part, as the SDK takes it. */
export type ReasoningPart = {
  type: 'reasoning';
  text: string;
  providerOptions: Record<string, Record<string, string>>;
};

function streamFailure(error: unknown): Error {
  if (error instanceof Error) return error;
  return new Error(typeof error === 'string' && error.length > 0 ? error : 'model stream failed');
}

function waitForAbort(signal: AbortSignal): Promise<never> {
  const aborted = new Promise<never>((_, reject) => {
    const abort = () => reject(signal.reason instanceof Error ? signal.reason : new Error('aborted'));
    if (signal.aborted) abort();
    else signal.addEventListener('abort', abort, { once: true });
  });
  void aborted.catch(() => undefined);
  return aborted;
}

export function streamTextDelta(part: StreamPart): string {
  if (typeof part.text === 'string' && part.text.length > 0) return part.text;
  return typeof part.delta === 'string' ? part.delta : '';
}

/** The call keeps the name the model gave unless it names a known tool, which maps back to its id. */
function streamToolCall(part: StreamPart, tools: readonly WireTool[]): ToolCallPart {
  const wire = part.toolName ?? '';
  const tool = tools.find((item) => item.wireName === wire || item.id === wire);
  return {
    kind: 'tool-call',
    requestId: part.toolCallId ?? wire,
    tool: { id: tool?.id ?? wire },
    arguments: (part.input ?? {}) as Json,
  };
}

function streamUsage(part: StreamPart): ProviderUsage | undefined {
  if (part.totalUsage === undefined) return undefined;
  const usage = providerUsageFromSdk(part.totalUsage);
  return Object.keys(usage).length === 0 ? undefined : usage;
}

/** The SDK's `inputTokens` already counts cache reads and writes, as `ProviderUsage` wants. */
export function providerUsageFromSdk(usage: LanguageModelUsage): ProviderUsage {
  return {
    ...(typeof usage.inputTokens === 'number' ? { inputTokens: usage.inputTokens } : {}),
    ...(typeof usage.outputTokens === 'number' ? { outputTokens: usage.outputTokens } : {}),
    ...(typeof usage.inputTokenDetails?.cacheReadTokens === 'number'
      ? { cacheReadInputTokens: usage.inputTokenDetails.cacheReadTokens }
      : {}),
    ...(typeof usage.inputTokenDetails?.cacheWriteTokens === 'number'
      ? { cacheWriteInputTokens: usage.inputTokenDetails.cacheWriteTokens }
      : {}),
    ...(typeof usage.outputTokenDetails?.reasoningTokens === 'number'
      ? { reasoningOutputTokens: usage.outputTokenDetails.reasoningTokens }
      : {}),
  };
}

export function mapFinishReason(reason: string | undefined): ModelFinish {
  switch (reason) {
    case 'stop':
      return 'complete';
    case 'length':
      return 'output-limit';
    case 'content-filter':
      return 'filtered';
    case 'tool-calls':
      return 'tool-calls';
    default:
      throw new Error(`Model stream finished unsuccessfully (${reason ?? 'missing'}).`);
  }
}

/** A failed call as the harness takes it. An HTTP error response means the provider refused before
 * generating (`effect: "none"`); rate limits, timeouts and 5xx are transient, other 4xx are not. A
 * failure with no status (network, mid-stream) is transient with an unknown effect. The wait a
 * response asks for (`retry-after-ms`, `retry-after`) is kept as `retryAfterMs`. API keys never
 * appear in the message. */
export function modelError(error: unknown, apiKey: string): ModelError {
  const message = error instanceof Error ? error.message : String(error);
  const safe = apiKey.length > 0 ? message.split(apiKey).join('[redacted]') : message;
  if (error instanceof ModelError) return new ModelError(safe, error);
  const source = error as {
    statusCode?: unknown;
    isRetryable?: unknown;
    responseBody?: unknown;
    responseHeaders?: unknown;
  } | null;
  const status = typeof source?.statusCode === 'number' ? source.statusCode : undefined;
  // What the provider said, when the SDK did not read it (a body other than its own error format,
  // e.g. `{"detail": "…"}`): kept in the message, without the key.
  const detail = providerDetail(source?.responseBody);
  if (detail && !safe.includes(detail)) {
    const told = apiKey.length > 0 ? detail.split(apiKey).join('[redacted]') : detail;
    return modelError(
      Object.assign(new Error(`${safe}: ${told.slice(0, 500)}`), {
        statusCode: status,
        isRetryable: source?.isRetryable,
        responseHeaders: source?.responseHeaders,
      }),
      apiKey,
    );
  }
  const transient =
    status !== undefined
      ? status === 408 || status === 409 || status === 429 || status >= 500
      : typeof source?.isRetryable === 'boolean'
        ? source.isRetryable
        : true;
  // The status stays readable after the harness keeps only the failure's message (e.g. "… (HTTP 404)").
  const said = status !== undefined && !safe.includes(`HTTP ${status}`) ? `${safe} (HTTP ${status})` : safe;
  const retryAfterMs = retryAfterOf(source?.responseHeaders);
  return new ModelError(said, {
    transient,
    effect: status !== undefined ? 'none' : 'unknown',
    ...(retryAfterMs !== undefined ? { retryAfterMs } : {}),
  });
}

/** The wait a response asks for, in ms: `retry-after-ms`, else `retry-after` in seconds or as an HTTP date. */
export function retryAfterOf(headers: unknown, now = Date.now()): number | undefined {
  if (typeof headers !== 'object' || headers === null) return undefined;
  const get = (name: string) => {
    const found = Object.entries(headers).find(([k]) => k.toLowerCase() === name)?.[1];
    return typeof found === 'string' && found.trim() !== '' ? found.trim() : undefined;
  };
  const ms = Number(get('retry-after-ms'));
  if (Number.isFinite(ms) && ms >= 0) return Math.round(ms);
  const after = get('retry-after');
  if (after === undefined) return undefined;
  const seconds = Number(after);
  if (Number.isFinite(seconds)) return seconds >= 0 ? Math.round(seconds * 1000) : undefined;
  const at = Date.parse(after);
  return Number.isNaN(at) ? undefined : Math.max(0, at - now);
}

/** A provider's message in an error body: `detail`, `error.message`, `message` or `error`. */
function providerDetail(body: unknown): string | undefined {
  if (typeof body !== 'string' || body.length === 0) return undefined;
  try {
    const v = JSON.parse(body) as { detail?: unknown; message?: unknown; error?: unknown };
    const e = v.error as { message?: unknown } | string | undefined;
    const found = v.detail ?? (typeof e === 'object' ? e?.message : e) ?? v.message;
    return typeof found === 'string' && found.length > 0 ? found : undefined;
  } catch {
    return undefined;
  }
}
