/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type AssistantMessage,
  canonicalJson,
  digest,
  type Json,
  type ModelAdapter,
  type ModelBindings,
  type ModelCapabilities,
  ModelError,
  type ModelFinish,
  type ModelFrame,
  type ProjectedRequest,
  type ProjectInput,
  type ProviderUsage,
  UnrepresentableInput,
} from '@molstar/ai-agent';
import type { Clock } from '@molstar/ai-harness';

/** One scripted model turn. */
export interface ModelStep {
  readonly text?: string;
  readonly reasoning?: string;
  readonly calls?: readonly { readonly tool: string; readonly arguments: Json; readonly requestId?: string }[];
  /** A full message instead of `text` and `calls`. */
  readonly message?: AssistantMessage;
  readonly finish?: ModelFinish;
  /** Omitted: estimated from the request and the message. `null`: the provider reports no usage. */
  readonly usage?: ProviderUsage | null;
  /** The attempt fails with a `ModelError`, after any frames. */
  readonly fail?: { readonly message?: string; readonly transient?: boolean; readonly retryAfterMs?: number };
  /** Waits this long (virtual time) before the result, after the deltas. */
  readonly waitMs?: number;
  /** Text deltas to stream (default: the text in one delta). */
  readonly deltas?: readonly string[];
}

export interface ScriptedModelOptions {
  readonly clock?: Clock;
  readonly capabilities?: Partial<ModelCapabilities>;
  /** Overrides the payload (to make a retry's request differ, for example). */
  readonly payload?: (input: ProjectInput, projection: number) => Json;
  /** Throws `UnrepresentableInput` when it returns a reason. */
  readonly unrepresentable?: (input: ProjectInput) => string | undefined;
}

export interface ScriptedModel {
  readonly bindings: ModelBindings;
  readonly adapter: ModelAdapter;
  /** Every projection, in order. */
  readonly projections: ProjectInput[];
  /** Every invocation's request, in order. */
  readonly invocations: ProjectedRequest[];
  /** Steps not yet played. */
  remaining(): number;
  /** Adds steps to the script. */
  push(...steps: ModelStep[]): void;
}

const DEFAULT_CAPABILITIES: ModelCapabilities = {
  contextWindowTokens: 200_000,
  maxOutputTokens: 8_000,
  defaultMaxOutputTokens: 4_000,
  outputLimitIncludesReasoning: false,
  appendsTools: true,
  cacheScope: true,
};

function messageOf(step: ModelStep, index: number): AssistantMessage {
  if (step.message !== undefined) return step.message;
  const parts: AssistantMessage['parts'][number][] = [];
  if (step.text !== undefined) parts.push({ kind: 'text', text: step.text });
  for (const [n, call] of (step.calls ?? []).entries())
    parts.push({
      kind: 'tool-call',
      requestId: call.requestId ?? `req-${index}-${n}`,
      tool: { id: call.tool },
      arguments: call.arguments,
    });
  return { role: 'assistant', parts };
}

/** A model adapter that plays a script, deterministic and instant in virtual time. */
export function scriptedModel(script: readonly ModelStep[] = [], options: ScriptedModelOptions = {}): ScriptedModel {
  const steps = [...script];
  const projections: ProjectInput[] = [];
  const invocations: ProjectedRequest[] = [];
  let played = 0;
  const capabilities = { ...DEFAULT_CAPABILITIES, ...options.capabilities };

  const adapter: ModelAdapter = {
    capabilitiesFor: async () => capabilities,
    async project(input) {
      const reason = options.unrepresentable?.(input);
      if (reason !== undefined) throw new UnrepresentableInput(reason);
      projections.push(input);
      const payload =
        options.payload?.(input, projections.length) ??
        ({
          // As an adapter sends them: a message's label is never sent.
          messages: input.messages.map(({ role, parts }) => ({ role, parts })) as unknown as Json,
          tools: input.tools.map((tool) => tool.tool.id),
          configuration: input.configuration as unknown as Json,
        } satisfies Json);
      return {
        payload,
        digest: await digest(payload),
        adapterRevision: 'scripted-1',
        inputTokens: Math.ceil(canonicalJson(payload).length / 4),
      };
    },
    async *invoke({ request, signal }): AsyncIterable<ModelFrame> {
      invocations.push(request);
      const step = steps.shift();
      const index = played++;
      if (step === undefined) throw new Error(`scriptedModel: no step left for invocation ${index + 1}`);
      if (step.reasoning !== undefined) yield { kind: 'reasoning-delta', text: step.reasoning };
      for (const text of step.deltas ?? (step.text !== undefined ? [step.text] : [])) {
        signal.throwIfAborted();
        yield { kind: 'text-delta', text };
      }
      if (step.waitMs !== undefined) {
        if (options.clock === undefined) throw new Error('scriptedModel: waitMs needs a clock');
        await options.clock.sleep(step.waitMs, signal);
      }
      signal.throwIfAborted();
      if (step.fail !== undefined)
        throw new ModelError(step.fail.message ?? 'scripted failure', {
          transient: step.fail.transient ?? false,
          effect: 'none',
          ...(step.fail.retryAfterMs !== undefined ? { retryAfterMs: step.fail.retryAfterMs } : {}),
        });
      const message = messageOf(step, index);
      const hasCalls = message.parts.some((part) => part.kind === 'tool-call');
      const usage: ProviderUsage =
        step.usage === null
          ? {}
          : (step.usage ?? {
              inputTokens: request.inputTokens,
              outputTokens: Math.ceil(canonicalJson(message as unknown as Json).length / 4),
            });
      yield { kind: 'result', message, finish: step.finish ?? (hasCalls ? 'tool-calls' : 'complete'), usage };
    },
  };

  return {
    adapter,
    bindings: { resolve: async () => adapter },
    projections,
    invocations,
    remaining: () => steps.length,
    push: (...more) => {
      steps.push(...more);
    },
  };
}
