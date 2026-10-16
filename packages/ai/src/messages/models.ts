/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AssistantMessage, Message, ToolDefinition } from './messages.js';
import type { Json } from './values.js';

/** Which model; the injected resolver owns routing and credentials. */
export interface ModelTarget {
  readonly provider: string;
  readonly model: string;
}

export type ReasoningParameters =
  | { readonly mode: 'off' }
  | { readonly mode: 'adaptive'; readonly effort?: string }
  | { readonly mode: 'budget'; readonly tokens: number };

/** Requests, not promises: a provider rejects what its model cannot do. */
export interface ModelParameters {
  readonly temperature?: number;
  readonly topP?: number;
  readonly maxOutputTokens?: number;
  readonly reasoning?: ReasoningParameters;
}

export interface ModelConfiguration {
  readonly target: ModelTarget;
  readonly parameters: ModelParameters;
  /** A prompt-cache hint; providers whose capabilities support it map it, others ignore it. */
  readonly cacheScopeId?: string;
}

export interface ModelCapabilities {
  readonly contextWindowTokens: number;
  /** An input cap for providers whose window does not include the output. */
  readonly maxInputTokens?: number;
  readonly maxOutputTokens?: number;
  /** The output limit used when `maxOutputTokens` is not set. */
  readonly defaultMaxOutputTokens?: number;
  readonly outputLimitIncludesReasoning: boolean;
  /** Whether a definition loaded later can be appended where it loads, not hoisted. */
  readonly appendsTools: boolean;
  /** Whether `cacheScopeId` is mapped to the provider's caching. */
  readonly cacheScope: boolean;
}

/** Provider-reported counts for one attempt. */
export interface ProviderUsage {
  /** Total input, cache reads and writes included. */
  readonly inputTokens?: number;
  readonly outputTokens?: number;
  readonly cacheReadInputTokens?: number;
  readonly cacheWriteInputTokens?: number;
  readonly reasoningOutputTokens?: number;
  readonly providerRequestId?: string;
}

export interface ProjectInput {
  readonly messages: readonly Message[];
  /** The available tools: enabled and loaded in the head. */
  readonly tools: readonly ToolDefinition[];
  readonly configuration: ModelConfiguration;
  readonly capabilities: ModelCapabilities;
}

/** A provider payload; the harness keeps its digest and never interprets it. */
export interface ProjectedRequest {
  readonly payload: Json;
  readonly digest: string;
  /** The provider's projection and replay rules: a retry refuses a request of another (harness.md §6). */
  readonly revision: string;
  /** An estimate of the payload's input tokens, for admission. */
  readonly inputTokens: number;
  /**
   * The payload's parts as the provider caches them (tools, system, then each message), a short
   * digest each: a request extends an earlier one, and can reuse its cache, when the earlier one's
   * segments are a prefix of its own. Where a provider gives them.
   */
  readonly segments?: readonly string[];
}

export type ModelFinish = 'complete' | 'tool-calls' | 'output-limit' | 'filtered';

export type ModelFrame =
  | { readonly kind: 'text-delta'; readonly text: string }
  | { readonly kind: 'reasoning-delta'; readonly text: string }
  /** Totals so far; the last one counts. */
  | { readonly kind: 'usage'; readonly usage: ProviderUsage }
  /** Exactly one, last. */
  | {
      readonly kind: 'result';
      readonly message: AssistantMessage;
      readonly finish: ModelFinish;
      readonly usage: ProviderUsage;
    };

/** A model provider (Anthropic, OpenAI, Google, an endpoint of one's own): its models' capabilities, a request
 * projected to its wire, and the call. */
export interface Provider {
  capabilities(configuration: ModelConfiguration): Promise<ModelCapabilities>;
  /** Pure: the same input gives the same payload and digest. */
  project(input: ProjectInput): Promise<ProjectedRequest>;
  invoke(input: { readonly request: ProjectedRequest; readonly signal: AbortSignal }): AsyncIterable<ModelFrame>;
}

/** The harness's way to its providers: one for a target. */
export interface ProviderBindings {
  resolve(target: ModelTarget): Promise<Provider>;
}

/** Input a provider cannot represent (a role where it may not stand, for example). */
export class UnrepresentableInput extends Error {
  readonly code = 'unrepresentable';
}

/**
 * A failed model call, as a provider reports it. By default not transient, with an unknown effect: a provider says
 * `effect: 'none'` when the request surely did nothing (so a retry is safe), and `retryAfterMs` when the provider said
 * how long to wait (`retry-after`).
 */
export class ModelError extends Error {
  readonly transient: boolean;
  readonly effect: 'none' | 'unknown';
  readonly retryAfterMs?: number;
  constructor(
    message: string,
    options: {
      readonly transient?: boolean;
      readonly effect?: 'none' | 'unknown';
      readonly retryAfterMs?: number;
    } = {},
  ) {
    super(message);
    this.name = 'ModelError';
    this.transient = options.transient ?? false;
    this.effect = options.effect ?? 'unknown';
    if (options.retryAfterMs !== undefined) this.retryAfterMs = options.retryAfterMs;
  }
}
