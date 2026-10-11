/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message, ModelConfiguration, ToolDefinition } from '@molstar/ai-agent';
import type { Limits, RunContext, RunResult, ToolHandler } from '@molstar/ai-harness';
import type { ToolRegistration } from './registry.js';

/**
 * A helper a tool runs under its own call (a nested run, harness.md §9): its model calls count in the turn's usage and
 * budget and show under the call in the trace. Its tools are its own (closures over the helper's state),
 * preauthorized, and seen by nothing else. `M` is the host's choice of model for it.
 */
export interface NestRequest<Api, M = unknown> {
  readonly messages: readonly Message[];
  readonly tools: readonly ToolRegistration<Api>[];
  /** The API its tools see. */
  readonly api: Api;
  readonly model: M;
  readonly limits: Limits;
}
export type Nest<Api, M = unknown> = (request: NestRequest<Api, M>) => Promise<RunResult>;

/** What a handler needs of the nested runs: to run one under the call. */
export interface NestedRuns<Api> {
  run(context: RunContext, request: NestRequest<Api, unknown>): Promise<RunResult>;
}

/** What the bindings need of the nested runs: to route a call to a helper's tool, and to run one. */
export interface NestedBinding<Api> extends NestedRuns<Api> {
  resolve(definition: ToolDefinition): ToolHandler | undefined;
}

/**
 * The tools of the nested runs in progress, by the call that runs each: bound once in the agent's bindings, each call
 * routed to its own run's handlers by the last call of its path.
 */
export class NestedTools<Api, M = unknown> implements NestedBinding<Api> {
  private readonly active = new Map<string, Map<string, ToolHandler>>();
  private readonly configure: (choice: M) => ModelConfiguration;

  constructor(configure: (choice: M) => ModelConfiguration) {
    this.configure = configure;
  }

  resolve(definition: ToolDefinition): ToolHandler | undefined {
    const id = definition.tool.id;
    if (![...this.active.values()].some((handlers) => handlers.has(id))) return undefined;
    return {
      invoke: (context) => {
        const handler = this.active.get(context.path.at(-1) ?? '')?.get(id);
        if (!handler) throw new Error(`${id} belongs to a helper that is not running.`);
        return handler.invoke(context);
      },
    };
  }

  async run(context: RunContext, request: NestRequest<Api, M>): Promise<RunResult> {
    const ids = request.tools.map((r) => r.definition.tool.id);
    this.active.set(context.callId, new Map(request.tools.map((r) => [r.definition.tool.id, r.handler(request.api)])));
    try {
      return await context.run({
        input: request.messages,
        tools: { enabled: ids, preauthorized: ids, definitions: request.tools.map((r) => r.definition) },
        model: this.configure(request.model),
        limits: request.limits,
      });
    } finally {
      this.active.delete(context.callId);
    }
  }
}
