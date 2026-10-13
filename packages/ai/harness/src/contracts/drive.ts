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
import type { HeadEdit } from './ports.js';
import type {
  AttemptFailure,
  Budget,
  CallId,
  ContextSize,
  Dimension,
  HarnessId,
  Head,
  HeadId,
  ItemKey,
  Limits,
  ScopeId,
  Seq,
  TurnConfiguration,
  TurnId,
} from './records.js';
import type { RunResult, RunSpec } from './runner.js';

/** Decisions and drives (harness.md §7). */

export interface ItemView {
  readonly key: ItemKey;
  readonly messages: readonly Message[];
  readonly options?: Json;
  readonly submittedAt: number;
  readonly interrupt: boolean;
}
export interface CallView {
  readonly callId: CallId;
  readonly toolId: string;
  readonly arguments: Json;
  readonly definition?: ToolDefinition;
}
export interface LastCompaction {
  readonly seq: Seq;
  readonly sizeAfter: ContextSize;
}

/**
 * Calls of the step that have not started, given a failed result with `text` instead of running:
 * `calls` names which (a drive's view of what depends on the failure); absent, all of them.
 */
export interface SkipAnswer {
  readonly text: string;
  readonly calls?: readonly CallId[];
}

/** Each kind: [payload, answer]. */
export interface Decisions {
  permission: [{ call: CallView }, 'allow' | { deny: string } | { park: { note?: Json } }];
  schedule: [{ calls: readonly CallView[] }, { stages: readonly (readonly CallId[])[]; maxConcurrency?: number }];
  budget: [{ dimension: Dimension; budget: Budget; need?: number }, { extend: Limits }];
  failure: [
    {
      call?: CallView;
      attempt?: number;
      failure: AttemptFailure;
      /** A tool's failure: the step's calls that have not started (what `skip` may name). */
      unstarted?: readonly CallView[];
    },
    (
      | { retry: { afterMs?: number } }
      | { replace: ModelConfiguration }
      | { report: { skip?: SkipAnswer } }
      | { compact: HeadEdit }
    ),
  ];
  'model-result': [
    { finish: 'output-limit' | 'filtered'; message: AssistantMessage },
    'accept' | { replace: ModelConfiguration },
  ];
  repair: [
    {
      call: CallView;
      reason: 'unknown-tool' | 'not-loaded' | 'not-enabled' | 'no-handler' | 'invalid-arguments';
      /** The step's other calls that have not started (what `skip` may name). */
      unstarted: readonly CallView[];
    },
    { report: { parts?: readonly ModelVisiblePart[]; contributions?: readonly Message[]; skip?: SkipAnswer } },
  ];
  boundary: [
    { size: ContextSize; completing: boolean; lastCompaction?: LastCompaction },
    'proceed' | { steer: readonly ItemKey[] } | { compact: HeadEdit } | 'end',
  ];
  interrupt: [{ item: ItemView; stepKind: 'model' | 'tools' }, 'abandon' | 'skip'];
  'next-turn': [
    { size: ContextSize; lastCompaction?: LastCompaction },
    { begin: { keys: readonly ItemKey[]; compact?: HeadEdit } } | 'release',
  ];
}
export type DecisionKind = keyof Decisions;
/** Legal for every kind. */
export type Exit = 'stop' | { fail: Json };
export interface Decision<K extends DecisionKind = DecisionKind> {
  readonly kind: K;
  readonly path: readonly CallId[];
  readonly payload: Decisions[K][0];
}
export type Answer<K extends DecisionKind> = Decisions[K][1] | Exit;

export interface DriveContext {
  readonly harnessId: HarnessId;
  readonly turnId?: TurnId;
  readonly depth: number;
  readonly budget: Budget;
  readonly head: Head;
  readonly configuration?: TurnConfiguration;
  /** At depth 0; a nested drive sees the inbox by passing up. */
  readonly inbox: readonly ItemView[];
  readonly signal: AbortSignal;
  readonly journal: { read(scope: ScopeId): Promise<{ head?: HeadId; messages: Message[] }> };
  /** A nested run inside the open step (how a drive calls a model to compact). */
  run(spec: RunSpec): Promise<RunResult>;
}

export interface Drive {
  decide<K extends DecisionKind>(
    decision: Decision<K>,
    context: DriveContext,
  ): Answer<K> | 'pass' | Promise<Answer<K> | 'pass'>;
}

export type DriveHandlers = {
  [K in DecisionKind]?: (
    decision: Decision<K>,
    context: DriveContext,
  ) => Answer<K> | 'pass' | Promise<Answer<K> | 'pass'>;
};
