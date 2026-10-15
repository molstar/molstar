/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, Message } from '@molstar/ai-agent';
import type { HarnessObservable } from './observe.js';
import type {
  ActiveTurn,
  BlobRef,
  CallRecord,
  Compaction,
  EntryId,
  Execution,
  HarnessEvent,
  HarnessId,
  Head,
  HeadId,
  HeadRecord,
  ItemKey,
  LaunchKey,
  Limits,
  PendingStep,
  Pointer,
  Progress,
  ScopeId,
  Segment,
  Seq,
  StepEvent,
  StepId,
  TurnConfiguration,
  TurnEnd,
  TurnId,
  TurnSummary,
} from './records.js';

/** The host's handles (harness.md §12, §15): apps and runners never share one. */

export interface HarnessHost {
  create(harnessId: HarnessId, input: { readonly key: string }): Promise<void>;
  remove(harnessId: HarnessId): Promise<void>;
  client(harnessId: HarnessId): HarnessClient;
  runnerPort(harnessId: HarnessId): RunnerPort;
}

export interface HarnessReader {
  readonly id: HarnessId;
  execution(): Promise<Execution>;
  heads(ids: readonly HeadId[]): Promise<Head[]>;
  pointers(scopes: readonly ScopeId[]): Promise<(Pointer | undefined)[]>;
  entries(ids: readonly EntryId[]): Promise<Message[]>;
  events(input: { readonly after?: Seq; readonly limit?: number }): Promise<HarnessEvent[]>;
  /** Every event of each step, in seq order (a parked step has several). */
  steps(ids: readonly StepId[]): Promise<StepEvent[]>;
  /** Ended turns, newest first, plus the active one. */
  turns(input: { readonly before?: Seq; readonly limit?: number }): Promise<TurnSummary[]>;
  blobs(refs: readonly BlobRef[]): Promise<Json[]>;
  /** The current seq, then each change; level-triggered. */
  watch(): HarnessObservable<Seq>;
}

export type StopInput = { readonly key: string; readonly reason?: Json } & (
  | { readonly turnId: TurnId }
  | { readonly item: ItemKey }
);

export interface HarnessClient extends HarnessReader {
  submit(input: {
    readonly key: ItemKey;
    readonly messages: readonly Message[];
    readonly options?: Json;
    readonly interrupt?: { readonly turnId: TurnId };
  }): Promise<{ seq: Seq; launch?: LaunchKey }>;
  withdraw(input: { readonly key: ItemKey }): Promise<{ seq: Seq; withdrawn: boolean }>;
  stop(input: StopInput): Promise<{ seq: Seq; launch?: LaunchKey; withdrawn?: boolean }>;
  resume(input: { readonly key: string; readonly turnId: TurnId }): Promise<{ seq: Seq; launch?: LaunchKey }>;
  reset(input: {
    readonly key: string;
    readonly pointers?: Readonly<Record<ScopeId, HeadId>>;
    readonly reason?: Json;
  }): Promise<{ seq: Seq; launch?: LaunchKey }>;
  check(): Promise<{ seq: Seq; launch?: LaunchKey }>;
}

export interface Lease {
  readonly harnessId: HarnessId;
  readonly runnerId: string;
  readonly until: number;
}
export interface Committed {
  readonly seq: Seq;
  readonly lease: Lease;
}

/** One pointer move, built by the runner with the shared head scheme (`buildMove`). */
export interface Move {
  readonly pointer: ScopeId;
  readonly from?: HeadId;
  readonly head: HeadRecord;
  /** Segments the head record names that are new. */
  readonly segments: readonly Segment[];
}
/** A head edit as a drive answers it: the new head, kept hashes and new messages, in order. */
export interface HeadEdit {
  readonly entries: readonly (EntryId | Message)[];
}

export interface BeginInput {
  readonly key: string;
  readonly turnId: TurnId;
  readonly step: StepId;
  readonly items: readonly ItemKey[];
  readonly configuration: TurnConfiguration;
  readonly messages: readonly Message[];
  readonly move: Move;
  readonly compaction?: readonly Compaction[];
  readonly calls: readonly CallRecord[];
  readonly next: PendingStep;
  readonly progress?: Progress;
}
export interface CommitBatch {
  readonly key: string;
  readonly turnId: TurnId;
  readonly step: StepId;
  readonly ending: 'settled' | 'parked' | 'interrupted' | 'abandoned';
  readonly calls: readonly CallRecord[];
  /** New messages; the host stores those whose hash is new. */
  readonly messages: readonly Message[];
  readonly moves: readonly Move[];
  readonly compaction?: readonly Compaction[];
  readonly consumed: readonly ItemKey[];
  readonly next?: PendingStep;
  readonly end?: TurnEnd;
  readonly park?: ActiveTurn['park'];
  /** The turn's limits after an `extend`. */
  readonly limits?: Limits;
  readonly progress?: Progress;
}

export interface RunnerPort extends HarnessReader {
  /** Head records and segments as stored, so a cold runner builds near-constant appends (harness.md §14). */
  headRecords(ids: readonly HeadId[]): Promise<HeadRecord[]>;
  segments(ids: readonly string[]): Promise<Segment[]>;
  /** The last compaction of a pointer, from the event log's index (the `boundary` payload, §7.1). */
  lastCompaction(pointer: ScopeId): Promise<{ seq: Seq; sizeAfter: Compaction['sizeAfter'] } | undefined>;
  claim(input: {
    readonly launch: LaunchKey;
    readonly runnerId: string;
    readonly schemaVersion: number;
  }): Promise<{ lease: Lease; execution: Execution }>;
  /** Before the commit that refers to them; a ref is `digest(body)`. */
  putBlobs(lease: Lease, blobs: readonly { readonly ref: BlobRef; readonly body: Json }[]): Promise<void>;
  begin(lease: Lease, input: BeginInput): Promise<Committed>;
  commit(lease: Lease, batch: CommitBatch): Promise<Committed>;
  /** Only when nothing else went out for `signalEveryMs`. */
  renew(lease: Lease, input?: { readonly progress?: Progress }): Promise<Lease>;
  release(lease: Lease, input: { readonly seenSeq: Seq }): Promise<{ seq: Seq; launch?: LaunchKey }>;
}

export type HarnessErrorCode =
  | 'invalid'
  | 'head-changed'
  | 'inbox-changed'
  | 'stale'
  | 'stale-schema'
  | 'pointer-busy'
  | 'unknown-outcome'
  | 'integrity';

/** The one operation error (harness.md §16). */
export class HarnessError extends Error {
  readonly code: HarnessErrorCode;
  readonly retry: 'same-key' | 'no';
  readonly detail?: Json;
  constructor(code: HarnessErrorCode, message: string, detail?: Json) {
    super(message);
    this.name = 'HarnessError';
    this.code = code;
    this.retry = code === 'unknown-outcome' ? 'same-key' : 'no';
    if (detail !== undefined) this.detail = detail;
  }
}
