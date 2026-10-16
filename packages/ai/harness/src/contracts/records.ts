/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, ModelConfiguration } from '@molstar/ai-agent';

/** The records of harness.md §2–§5: what the harness holds in memory and the event log records. */

export type HarnessId = string;
export type ScopeId = string;
export type HeadId = string;
/** `hashMessage(message)`. */
export type EntryId = string;
export type TurnId = string;
export type StepId = string;
export type CallId = string;
export type AttemptId = string;
/** `digest(body)`. */
export type BlobRef = string;
export type Seq = number;

/** The top-level pointer's scope. */
export const MAIN: ScopeId = 'main';

/** A head: its entries' hashes in order (a hash may repeat); `id` digests them, so a head is content-addressed. */
export interface Head {
  readonly id: HeadId;
  readonly entries: readonly EntryId[];
}
export interface Pointer {
  readonly scope: ScopeId;
  readonly head: HeadId;
}

export interface ToolSet {
  readonly enabled: readonly string[];
  readonly preauthorized: readonly string[];
}
export interface TurnConfiguration {
  readonly model: ModelConfiguration;
  readonly tools: ToolSet;
  readonly limits: Limits;
}

export interface Usage {
  readonly modelCalls: number;
  readonly toolCalls: number;
  readonly uncachedInput: number;
  readonly cacheRead: number;
  readonly cacheWrite: number;
  readonly output: number;
  readonly reasoning: number;
  readonly cost?: number;
  readonly estimated?: boolean;
}
export type Dimension =
  | 'modelCalls'
  | 'toolCalls'
  | 'uncachedInput'
  | 'cacheRead'
  | 'cacheWrite'
  | 'output'
  | 'reasoning'
  | 'cost'
  | 'endBy';
export type Limits = Partial<Record<Dimension, number>>;
export interface Budget {
  readonly limits: Limits;
  readonly consumed: Usage;
}

export interface PlannedCall {
  readonly callId: CallId;
  readonly key: string;
  readonly requestId: string;
  readonly toolId: string;
  readonly arguments: Json;
  /** The entry holding the definition the model saw; absent: unbound. */
  readonly definition?: EntryId;
  readonly unbound?: 'unknown-tool' | 'not-loaded' | 'not-enabled' | 'no-handler' | 'invalid-arguments';
}
export type PendingStep =
  | { readonly id: StepId; readonly kind: 'model' }
  | { readonly id: StepId; readonly kind: 'tools'; readonly message: BlobRef; readonly plan: readonly PlannedCall[] };

export interface ActiveTurn {
  readonly id: TurnId;
  readonly configuration: TurnConfiguration;
  readonly budget: Budget;
  readonly pending: PendingStep;
  readonly park?: { readonly path: readonly CallId[]; readonly note?: Json };
  /** Running totals, nested runs included. */
  readonly usage: Usage;
}

/** The reader's snapshot of the harness (harness.md §4): the log's last seq, `main`, and the active turn. */
export interface Execution {
  readonly harnessId: HarnessId;
  readonly seq: Seq;
  readonly main: HeadId;
  readonly turn?: ActiveTurn;
}

export interface AttemptFailure {
  readonly code:
    | 'provider'
    | 'timeout'
    | 'stream-protocol'
    | 'invalid-contribution'
    | 'capacity'
    | 'unrepresentable'
    | 'aborted'
    | 'handler'
    | 'integrity';
  readonly message: string;
  readonly transient: boolean;
  readonly effect: 'none' | 'applied' | 'unknown';
  /** How long the provider asked to wait before another attempt (a model's `retry-after`), if it said. */
  readonly retryAfterMs?: number;
}
export interface AttemptRecord {
  readonly id: AttemptId;
  readonly startedAt: number;
  readonly endedAt: number;
  readonly outcome: 'succeeded' | 'failed' | 'aborted';
  readonly failure?: AttemptFailure;
  readonly usage: Usage;
}
export interface CallRecord {
  readonly id: CallId;
  readonly key: string;
  readonly parentCallId?: CallId;
  readonly kind: 'model' | 'tool';
  readonly request: BlobRef;
  readonly attempts: readonly AttemptRecord[];
  readonly result?: BlobRef;
}

export interface ContextSize {
  readonly inputTokens: number;
  readonly reservedOutput: number;
  readonly window: number;
}
export interface TurnEnd {
  readonly status: 'completed' | 'failed' | 'stopped' | 'reset';
  readonly reason?: Json;
}
export interface TurnSummary {
  readonly id: TurnId;
  readonly status: 'active' | TurnEnd['status'];
  /** The turn's begin step. */
  readonly firstStep: StepId;
  readonly lastStep?: StepId;
  readonly usage: Usage;
  readonly endedAt?: number;
}

export interface Compaction {
  readonly pointer: ScopeId;
  readonly sizeAfter: ContextSize;
}
export interface StepEvent {
  readonly kind: 'step';
  readonly seq: Seq;
  readonly at: number;
  readonly turnId: TurnId;
  readonly step: StepId;
  readonly stepKind: 'begin' | 'model' | 'tools';
  readonly ending: 'settled' | 'parked' | 'abandoned';
  readonly calls: readonly CallRecord[];
  readonly moves: readonly { readonly pointer: ScopeId; readonly from?: HeadId; readonly to: HeadId }[];
  readonly compaction?: readonly Compaction[];
  readonly next?: PendingStep;
  readonly end?: TurnEnd;
  readonly summary?: TurnSummary;
  /** On a begin step: the turn's configuration. */
  readonly begun?: { readonly configuration: TurnConfiguration };
  readonly limits?: Limits;
}
/** What happened outside a step: a parked turn stopped or resumed, a reset (harness.md §5). */
export interface ControlEvent {
  readonly kind: 'control';
  readonly seq: Seq;
  readonly at: number;
  readonly control: 'stopped' | 'resumed' | 'reset';
  readonly turnId?: TurnId;
  readonly summary?: TurnSummary;
  readonly detail: Json;
}
export type HarnessEvent = StepEvent | ControlEvent;
