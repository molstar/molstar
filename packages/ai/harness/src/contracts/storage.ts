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
  BlobRef,
  Execution,
  HarnessEvent,
  HarnessId,
  HeadRecord,
  Pointer,
  ScopeId,
  Segment,
  Seq,
  StepEvent,
  StepId,
  TurnId,
} from './records.js';

/**
 * The storage port (harness.md §14): what a backend implements. The protocol logic is the harness's; a
 * backend stores and indexes records. Everything but the Execution and pointers is written once: a
 * `put` that would change a stored message, head record, segment or event fails `integrity`.
 */
export interface HarnessStorage {
  /** Serializable per harness: the callback's writes commit together, or none do. */
  transact<T>(harnessId: HarnessId, run: (tx: StorageTx) => Promise<T>): Promise<T>;
  /** Outside transactions, before the commit that refers to them. */
  putBlobs(harnessId: HarnessId, blobs: readonly { readonly ref: BlobRef; readonly body: Json }[]): Promise<void>;
  getBlobs(harnessId: HarnessId, refs: readonly BlobRef[]): Promise<(Json | undefined)[]>;
  /** The Execution's seq now and on each change. */
  watch(harnessId: HarnessId): HarnessObservable<Seq>;
  /** Drops a harness's records (after its tombstone is written). */
  drop(harnessId: HarnessId): Promise<void>;
}

export interface RecordStore<T> {
  get(ids: readonly string[]): Promise<(T | undefined)[]>;
  put(id: string, value: T): void;
}

export interface EventLog {
  append(event: HarnessEvent): void;
  after(seq: Seq, limit: number): Promise<HarnessEvent[]>;
  byStep(ids: readonly StepId[]): Promise<StepEvent[]>;
  byTurn(turnId: TurnId, after: Seq | undefined, limit: number): Promise<HarnessEvent[]>;
  /** Turn-ending events, newest first, before `before`. */
  endedTurns(before: Seq | undefined, limit: number): Promise<HarnessEvent[]>;
  byKey(key: string): Promise<HarnessEvent[]>;
  lastCompaction(pointer: ScopeId): Promise<StepEvent | undefined>;
}

/** Store-level records: tombstones and create keys live outside any harness. */
export interface StoreTx {
  tombstone(harnessId: HarnessId): Promise<boolean>;
  setTombstone(harnessId: HarnessId): void;
}

export interface StorageTx {
  readonly execution: RecordStore<Execution>;
  readonly events: EventLog;
  readonly messages: RecordStore<Message>;
  readonly pointers: RecordStore<Pointer>;
  readonly heads: RecordStore<HeadRecord>;
  readonly segments: RecordStore<Segment>;
  /** Tombstones (harness.md §12, create and remove). */
  readonly store: StoreTx;
}
