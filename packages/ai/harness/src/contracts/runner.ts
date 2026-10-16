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
  Message,
  ModelConfiguration,
  ModelVisiblePart,
  ToolDefinition,
} from '@molstar/ai-agent';
import type { Drive } from './drive.js';
import type { HeadEdit } from './harness.js';
import type { CallId, HarnessId, HeadId, Limits, ScopeId, ToolSet, TurnId, Usage } from './records.js';

/** What the loop runs with (harness.md §6–§9): time and ids, injected; tool handlers and their contexts. */

/** Time, injected so every retry and deadline runs on a simulated clock in tests. */
export interface Clock {
  now(): number;
  /** Resolves after `ms`, or rejects with the signal's reason when it aborts. */
  sleep(ms: number, signal?: AbortSignal): Promise<void>;
}

/** Fresh ids for turns, steps, calls and attempts. */
export interface IdSource {
  next(kind: 'turn' | 'step' | 'call' | 'attempt'): string;
}

export interface ToolBindings {
  resolve(definition: ToolDefinition): ToolHandler | undefined;
}

export type ToolFrame =
  | { readonly kind: 'progress'; readonly value: Json }
  | { readonly kind: 'usage'; readonly usage: Partial<Usage> }
  | {
      readonly kind: 'result';
      readonly status: 'succeeded' | 'failed';
      readonly parts: readonly ModelVisiblePart[];
      readonly contributions?: readonly Message[];
      readonly effect?: 'none' | 'applied' | 'unknown';
      readonly transient?: boolean;
    };

export interface ToolHandler {
  invoke(context: RunContext): AsyncIterable<ToolFrame>;
}

/** A tool's own context (harness.md §8): reads see this step's staged moves; writes are staged. */
export interface ToolJournal {
  read(scope: ScopeId): Promise<{ head?: HeadId; messages: Message[] }>;
  /** Throws `pointer-busy` for another writer's pointer and for `main`; a failed write fails the attempt. */
  append(scope: ScopeId, messages: readonly Message[]): void;
  edit(scope: ScopeId, edit: HeadEdit): void;
}

export interface RunSpec {
  readonly input: readonly Message[];
  readonly tools: ToolSet & { readonly definitions: readonly ToolDefinition[] };
  readonly model: ModelConfiguration;
  /** Where the context starts: a pointer (continued and moved), a head (a fork), or neither (empty). */
  readonly pointer?: ScopeId;
  readonly from?: HeadId;
  readonly limits: Limits;
  readonly drive?: Drive;
}
export interface RunResult {
  readonly status: 'completed' | 'failed';
  readonly output?: AssistantMessage;
  readonly usage: Usage;
  readonly failure?: Json;
}

export interface RunContext {
  readonly harnessId: HarnessId;
  readonly turnId: TurnId;
  readonly callId: CallId;
  readonly path: readonly CallId[];
  readonly key: string;
  /** `possible` on a retry: an earlier attempt may have had its effect. */
  readonly prior: 'none' | 'possible';
  readonly attempt: number;
  readonly definition: ToolDefinition;
  readonly arguments: Json;
  /** This call's own signal. */
  readonly signal: AbortSignal;
  readonly journal: ToolJournal;
  run(spec: RunSpec): Promise<RunResult>;
}
