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
  ModelFinish,
  ModelFrame,
  ProviderUsage,
  ToolCallPart,
} from '../../messages/index.js';
import { ModelError } from '../../messages/index.js';
import { jsonSchema, type LanguageModelUsage, streamText, type ToolSet, UnsupportedFunctionalityError } from 'ai';
import type { Payload } from './provider.js';

// The one stream mapper (`_spec/models.md`, "The providers"): the retained payload goes to `streamText` as it is (its tools' schemas
// wrapped by `jsonSchema`), the stream's parts become frames, and failures `ModelError`s.

/** What a provider keeps of a response for replay (signatures, encrypted reasoning, item ids): on the calls. */
export interface ReplayState {
  reasoning(part: StreamPart): void;
  toolCall(call: ToolCallPart, part: StreamPart, index: number): ToolCallPart;
  finish(calls: readonly ToolCallPart[]): readonly ToolCallPart[];
}

export const noReplay = (): ReplayState => ({ reasoning() {}, toolCall: (call) => call, finish: (calls) => calls });

type Model = Parameters<typeof streamText>[0]['model'];

/** The payload's call arguments: the identity, the model id and the tool ids are the provider's, not the SDK's. */
function callArguments(payload: Payload, model: Model, signal: AbortSignal, onError: (error: unknown) => void) {
  const { provider: _p, revision: _r, model: _m, toolIds: _i, appended: _a, tools, ...rest } = payload;
  const sdkTools: ToolSet = {};
  for (const [name, tool] of Object.entries(tools ?? {})) {
    sdkTools[name] = {
      description: tool.description,
      inputSchema: jsonSchema(tool.inputSchema as Parameters<typeof jsonSchema>[0]),
      ...(tool.providerOptions === undefined ? {} : { providerOptions: tool.providerOptions }),
    };
  }
  return {
    ...(rest as unknown as Omit<Parameters<typeof streamText>[0], 'model'>),
    ...(tools === undefined ? {} : { tools: sdkTools }),
    model,
    abortSignal: signal,
    maxRetries: 0,
    onError({ error }: { error: unknown }) {
      onError(error);
    },
  } as Parameters<typeof streamText>[0];
}

export async function* streamPayload(input: {
  readonly payload: Payload;
  readonly signal: AbortSignal;
  readonly languageModel: Model;
  readonly replay: ReplayState;
}): AsyncIterable<ModelFrame> {
  const { payload, replay } = input;
  let streamError: unknown;
  const result = streamText(
    callArguments(payload, input.languageModel, input.signal, (error) => {
      streamError = error;
    }),
  );
  const text: string[] = [];
  const toolCalls: ToolCallPart[] = [];
  let finish: ModelFinish | undefined;
  let usage: ProviderUsage = {};
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
        toolCalls.push(replay.toolCall(streamToolCall(typed, payload.toolIds), typed, toolCalls.length));
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
function streamToolCall(part: StreamPart, toolIds: Readonly<Record<string, string>>): ToolCallPart {
  const wire = part.toolName ?? '';
  return {
    kind: 'tool-call',
    requestId: part.toolCallId ?? wire,
    tool: { id: Object.hasOwn(toolIds, wire) ? (toolIds[wire] as string) : wire },
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
  // The SDK refuses what the provider cannot take before any request: retrying cannot help.
  if (UnsupportedFunctionalityError.isInstance(error))
    return new ModelError(safe, { transient: false, effect: 'none' });
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
