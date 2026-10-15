/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, Message, ModelVisiblePart, ToolDefinition } from '@molstar/ai-agent';
import type { CallId, RunContext, ToolBindings, ToolFrame, ToolHandler } from '@molstar/ai-harness';
import { toolDefinition } from './messages.js';

/** What a scripted tool returns; becomes the one `result` frame. */
export interface ToolOutcome {
  readonly status?: 'succeeded' | 'failed';
  readonly text?: string;
  readonly parts?: readonly ModelVisiblePart[];
  readonly contributions?: readonly Message[];
  readonly transient?: boolean;
  readonly effect?: 'none' | 'applied' | 'unknown';
}

/** A scripted tool: an outcome (sync or async), or a raw frame stream for protocol cases. */
export type ToolRun = (context: RunContext) => ToolOutcome | Promise<ToolOutcome> | AsyncIterable<ToolFrame>;

export interface ToolSpec {
  readonly definition?: ToolDefinition;
  /** Omitted: the definition binds to no handler (`no-handler`). */
  readonly run?: ToolRun;
}

export interface ToolCallRecord {
  readonly toolId: string;
  readonly callId: CallId;
  readonly key: string;
  readonly prior: 'none' | 'possible';
  readonly attempt: number;
  readonly arguments: Json;
  readonly path: readonly CallId[];
  readonly signal: AbortSignal;
  ended: boolean;
}

export interface ScriptedTools {
  readonly bindings: ToolBindings;
  readonly definitions: readonly ToolDefinition[];
  definition(id: string): ToolDefinition;
  /** Every invocation, in start order. */
  readonly calls: ToolCallRecord[];
  callsOf(id: string): ToolCallRecord[];
}

function isStream(value: unknown): value is AsyncIterable<ToolFrame> {
  return typeof value === 'object' && value !== null && Symbol.asyncIterator in value;
}

/** Tool bindings over scripted handlers; an echo tool when `run` is `true`. */
export function scriptedTools(specs: Readonly<Record<string, ToolSpec | ToolRun>>): ScriptedTools {
  const calls: ToolCallRecord[] = [];
  const definitions: ToolDefinition[] = [];
  const handlers = new Map<string, ToolHandler>();
  for (const [id, raw] of Object.entries(specs)) {
    const spec: ToolSpec = typeof raw === 'function' ? { run: raw } : raw;
    definitions.push(spec.definition ?? toolDefinition(id));
    const run = spec.run;
    if (run === undefined) continue;
    handlers.set(id, {
      async *invoke(context) {
        const record: ToolCallRecord = {
          toolId: id,
          callId: context.callId,
          key: context.key,
          prior: context.prior,
          attempt: context.attempt,
          arguments: context.arguments,
          path: context.path,
          signal: context.signal,
          ended: false,
        };
        calls.push(record);
        try {
          const outcome = run(context);
          if (isStream(outcome)) {
            yield* outcome;
            return;
          }
          const result = await outcome;
          yield {
            kind: 'result',
            status: result.status ?? 'succeeded',
            parts: result.parts ?? [{ kind: 'text', text: result.text ?? `${id} done` }],
            ...(result.contributions !== undefined ? { contributions: result.contributions } : {}),
            ...(result.transient !== undefined ? { transient: result.transient } : {}),
            ...(result.effect !== undefined ? { effect: result.effect } : {}),
          };
        } finally {
          record.ended = true;
        }
      },
    });
  }
  return {
    definitions,
    calls,
    definition(id) {
      const found = definitions.find((definition) => definition.tool.id === id);
      if (found === undefined) throw new Error(`scriptedTools: no tool ${id}`);
      return found;
    },
    callsOf: (id) => calls.filter((call) => call.toolId === id),
    bindings: { resolve: (definition) => handlers.get(definition.tool.id) },
  };
}
