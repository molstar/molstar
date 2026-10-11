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
  ModelBindings,
  ModelConfiguration,
  ModelVisiblePart,
  ToolDefinition,
} from '@molstar/ai-agent';
import type { Drive } from './drive.js';
import type { HarnessObservable, Observation } from './observe.js';
import type { HeadEdit, RunnerPort } from './ports.js';
import type {
  CallId,
  HarnessId,
  HeadId,
  Item,
  LaunchKey,
  Limits,
  ScopeId,
  Seq,
  ToolSet,
  TurnConfiguration,
  TurnId,
  Usage,
} from './records.js';

/** Time, injected so every lease, retry and deadline runs on a simulated clock in tests. */
export interface Clock {
  now(): number;
  /** Resolves after `ms`, or rejects with the signal's reason when it aborts. */
  sleep(ms: number, signal?: AbortSignal): Promise<void>;
}

/** Fresh ids for turns, steps, calls, attempts and runners; collision-resistant across workers. */
export interface IdSource {
  next(kind: 'turn' | 'step' | 'call' | 'attempt' | 'runner' | 'launch' | 'key'): string;
}

export interface RunnerOptions {
  readonly port: RunnerPort;
  readonly runnerId: string;
  readonly models: ModelBindings;
  readonly tools: ToolBindings;
  readonly drive: Drive;
  configure(items: readonly Item[]): TurnConfiguration | Promise<TurnConfiguration>;
  price?(configuration: ModelConfiguration, usage: Usage): number;
  readonly clock?: Clock;
  readonly ids?: IdSource;
  /** Default: `SCHEMA_VERSION`; another value only in tests (a runner of another version). */
  readonly schemaVersion?: number;
  readonly limits?: {
    readonly signalEveryMs?: number;
    readonly handOffAfterMs?: number;
    readonly depth?: number;
    readonly cacheBytes?: number;
  };
}

export interface RunnerExit {
  readonly reason: 'released' | 'handed-off' | 'lost' | 'stale-launch' | 'stale-schema' | 'aborted';
  readonly seq?: Seq;
}

export interface Runner {
  /** Claims with the launch key, serves turns, then releases or hands off. */
  run(launch: LaunchKey): Promise<RunnerExit>;
  observe(): HarnessObservable<Observation>;
  /** Local only: stop working now. Not a durable Stop (that is the client's `stop`). */
  abort(reason?: string): void;
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
  /** Throws `pointer-busy` for another writer's pointer and for `main`. */
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
  readonly prior: 'none' | 'possible';
  readonly attempt: number;
  readonly definition: ToolDefinition;
  readonly arguments: Json;
  /** This call's own signal. */
  readonly signal: AbortSignal;
  readonly journal: ToolJournal;
  run(spec: RunSpec): Promise<RunResult>;
}
