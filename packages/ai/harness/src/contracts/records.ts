/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, Message, ModelConfiguration } from '@molstar/ai-agent';

/** The records of harness.md §2–§4 and api.md: what the host stores and the event log holds. */

export type HarnessId = string;
export type ScopeId = string;
export type HeadId = string;
/** `hashMessage(message)`. */
export type EntryId = string;
export type TurnId = string;
export type StepId = string;
export type CallId = string;
export type AttemptId = string;
export type ItemKey = string;
export type LaunchKey = string;
/** `digest(body)`. */
export type BlobRef = string;
export type Seq = number;

/** The top-level pointer's scope. */
export const MAIN: ScopeId = 'main';

/** The schema version this harness writes; a turn written by another one is blocked (harness.md §4). */
export const SCHEMA_VERSION = 1;

/** A head as read: its entries' hashes in order (a hash may repeat). */
export interface Head {
  readonly id: HeadId;
  readonly entries: readonly EntryId[];
}
/** How a head is stored: its segments, content-addressed (`id` digests the segment ids). */
export interface HeadRecord {
  readonly id: HeadId;
  readonly segments: readonly string[];
  readonly length: number;
}
/** An immutable chunk of entry hashes, content-addressed. */
export interface Segment {
  readonly id: string;
  readonly entries: readonly EntryId[];
}
export interface Pointer {
  readonly scope: ScopeId;
  readonly head: HeadId;
}

export interface Item {
  readonly key: ItemKey;
  readonly messages: readonly Message[];
  readonly options?: Json;
  readonly submittedAt: number;
  readonly interrupt?: { readonly turnId: TurnId };
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
  readonly schemaVersion: number;
  readonly configuration: TurnConfiguration;
  readonly budget: Budget;
  readonly pending: PendingStep;
  readonly park?: { readonly path: readonly CallId[]; readonly note?: Json };
  readonly stopRequested?: { readonly reason?: Json; readonly at: number };
  readonly blocked?: 'stale-schema';
  /** Running totals, nested runs included. */
  readonly usage: Usage;
}

export interface Progress {
  readonly step?: StepId;
  readonly phase: 'model' | 'tools' | 'decision' | 'idle';
  readonly since: number;
  readonly detail?: Json;
}
export type RunnerState =
  | { readonly phase: 'none' }
  | { readonly phase: 'launching'; readonly key: LaunchKey; readonly until: number }
  | {
      readonly phase: 'live';
      readonly runnerId: string;
      readonly leaseUntil: number;
      readonly signalAt: number;
      readonly progress?: Progress;
    };

export interface Execution {
  readonly harnessId: HarnessId;
  readonly seq: Seq;
  readonly schemaVersion: number;
  readonly main: HeadId;
  readonly inbox: readonly Item[];
  readonly runner: RunnerState;
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
  readonly startedBy: readonly ItemKey[];
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
  readonly ending: 'settled' | 'parked' | 'interrupted' | 'abandoned';
  readonly calls: readonly CallRecord[];
  readonly moves: readonly { readonly pointer: ScopeId; readonly from?: HeadId; readonly to: HeadId }[];
  readonly compaction?: readonly Compaction[];
  readonly consumed: readonly ItemKey[];
  readonly next?: PendingStep;
  readonly end?: TurnEnd;
  readonly summary?: TurnSummary;
  readonly begun?: { readonly items: readonly ItemKey[]; readonly configuration: TurnConfiguration };
  readonly limits?: Limits;
  readonly commitKey: string;
  /** `digest` of the commit's input, so a reused key with other input is refused. */
  readonly inputDigest?: string;
}
export interface ControlEvent {
  readonly kind: 'control';
  readonly seq: Seq;
  readonly at: number;
  readonly control:
    | 'submitted'
    | 'withdrawn'
    | 'stop-requested'
    | 'stopped'
    | 'resumed'
    | 'reset'
    | 'claimed'
    | 'blocked'
    | 'released'
    | 'launched';
  readonly key?: string;
  /** `digest` of the operation's input, so a reused key with other input is refused. */
  readonly inputDigest?: string;
  readonly turnId?: TurnId;
  readonly summary?: TurnSummary;
  readonly detail: Json;
}
export type HarnessEvent = StepEvent | ControlEvent;
