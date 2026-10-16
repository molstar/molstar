/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, Message, ModelBindings, ModelConfiguration } from '@molstar/ai-agent';
import type { Drive } from './drive.js';
import type { HarnessObservable, Observation } from './observe.js';
import type {
  ActiveTurn,
  BlobRef,
  EntryId,
  Execution,
  HarnessEvent,
  HarnessId,
  Head,
  HeadId,
  Pointer,
  ScopeId,
  Seq,
  StepEvent,
  StepId,
  TurnConfiguration,
  TurnEnd,
  TurnId,
  TurnSummary,
  Usage,
} from './records.js';
import type { Clock, IdSource, ToolBindings } from './runner.js';

/** The harness (harness.md §4–§6): one in-memory loop, its reads, and its observations. */

/** What the trace and the Inspector read (harness.md §3, §13). Reads are async so a reader may be remote. */
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

/** How a turn left the loop: it ended, or it parked (`resume` goes on with it). */
export type TurnOutcome =
  | { readonly turnId: TurnId; readonly end: TurnEnd }
  | { readonly turnId: TurnId; readonly parked: NonNullable<ActiveTurn['park']> };

/** A turn the loop runs: its id (known at once, before its first step), and its outcome. */
export interface TurnHandle {
  readonly turnId: TurnId;
  /** Settles when the turn ends or parks; it does not reject for a failed turn (that is an outcome). */
  readonly done: Promise<TurnOutcome>;
}

export interface HarnessOptions {
  /** The reader's id (default: a fresh one). */
  readonly id?: HarnessId;
  readonly models: ModelBindings;
  readonly tools: ToolBindings;
  readonly drive: Drive;
  /** A turn's configuration when `turn` is given none. */
  configure?(messages: readonly Message[]): TurnConfiguration;
  price?(configuration: ModelConfiguration, usage: Usage): number;
  readonly clock?: Clock;
  readonly ids?: IdSource;
  readonly limits?: { readonly depth?: number };
}

export interface Harness extends HarnessReader {
  /** Begins a turn with these messages; refused (`invalid`) while a turn is active, running or parked. */
  turn(input: { readonly messages: readonly Message[]; readonly configuration?: TurnConfiguration }): TurnHandle;
  /**
   * Ends the active turn as stopped: a running step is aborted and recorded abandoned, with its attempts; a parked
   * turn ends at once. Nothing happens without an active turn.
   */
  stop(reason?: Json): void;
  /** Re-runs a parked turn's step, asking the drive again; refused (`invalid`) unless a turn is parked. */
  resume(): TurnHandle;
  /**
   * Ends any turn as reset (a running step's writes are dropped) and moves pointers: the given heads, else `main` to
   * its head at the end of the last completed turn. Deletes nothing.
   */
  reset(heads?: Readonly<Record<ScopeId, HeadId>>, reason?: Json): void;
  /** Live observations (deltas, attempts, decisions) and each event as it is written (harness.md §13). */
  observe(): HarnessObservable<Observation>;
}

/** A head edit as a drive answers it: the new head, kept hashes and new messages, in order. */
export interface HeadEdit {
  readonly entries: readonly (EntryId | Message)[];
}

export type HarnessErrorCode = 'invalid' | 'pointer-busy';

/** The one operation error (harness.md §16): a call the harness refuses. */
export class HarnessError extends Error {
  readonly code: HarnessErrorCode;
  readonly detail?: Json;
  constructor(code: HarnessErrorCode, message: string, detail?: Json) {
    super(message);
    this.name = 'HarnessError';
    this.code = code;
    if (detail !== undefined) this.detail = detail;
  }
}
