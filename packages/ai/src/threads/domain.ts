/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ThreadCommandId, ThreadEntryId, ThreadId, ThreadOwnerId } from './ids.js';
import type { ThreadLimits } from './limits.js';

export const THREAD_PAYLOAD_VERSION = 1;

export type ThreadPayloadEnvelope = Readonly<{
  readonly codec: string;
  readonly version: typeof THREAD_PAYLOAD_VERSION;
  readonly value: unknown;
}>;

export type ThreadAbortReason = 'caller-interrupted' | 'inactivity' | 'limit-exceeded' | 'host-closed';

export type ThreadEntryState =
  | Readonly<{
      readonly kind: 'streaming';
      readonly nextOffset: number;
      readonly lastActivityAt: number;
    }>
  | Readonly<{
      readonly kind: 'settled';
      readonly payload: ThreadPayloadEnvelope;
      readonly settledAt: number;
    }>
  | Readonly<{
      readonly kind: 'aborted';
      readonly reason: ThreadAbortReason;
      readonly endedAt: number;
    }>;

export type ThreadEntry = Readonly<{
  readonly entryId: ThreadEntryId;
  readonly threadId: ThreadId;
  readonly position: number;
  readonly descriptor: ThreadPayloadEnvelope;
  readonly state: ThreadEntryState;
}>;

export type ThreadDeltaBatch = Readonly<{
  readonly entryId: ThreadEntryId;
  readonly start: number;
  readonly endExclusive: number;
  readonly payload: ThreadPayloadEnvelope;
}>;

export type ThreadSnapshot = Readonly<{
  readonly threadId: ThreadId;
  readonly ownerId: ThreadOwnerId;
  readonly nextPosition: number;
  readonly lifecycle: 'active' | 'closing' | 'removed';
}>;

export type ThreadWriteResult = Readonly<{
  readonly thread: ThreadSnapshot;
  readonly entry: ThreadEntry;
}>;

export type ThreadDeltaResult = Readonly<{
  readonly thread: ThreadSnapshot;
  readonly entry: ThreadEntry;
  readonly nextOffset: number;
}>;

export type ThreadAppendDeltaResult =
  | Readonly<{
      readonly kind: 'appended';
      readonly thread: ThreadSnapshot;
      readonly entry: ThreadEntry;
      readonly nextOffset: number;
    }>
  | Readonly<{
      readonly kind: 'aborted';
      readonly thread: ThreadSnapshot;
      readonly entry: ThreadEntry;
      readonly reason: Extract<ThreadAbortReason, 'limit-exceeded'>;
    }>;

export type ThreadWindow = Readonly<{
  readonly threadId: ThreadId;
  readonly entries: readonly ThreadEntry[];
  readonly start: number;
  readonly endExclusive: number;
  readonly hasMore: boolean;
}>;

export type ThreadDeltaPage = Readonly<{
  readonly entryId: ThreadEntryId;
  readonly batches: readonly ThreadDeltaBatch[];
  readonly nextOffset: number;
  readonly hasMore: boolean;
  readonly gap?: Readonly<{ readonly afterOffset: number }>;
}>;

export type ThreadFailure = Readonly<{
  readonly status: 'conflict' | 'invalid' | 'not-found' | 'closed';
  readonly reason: string;
}>;

export type ThreadCommandResult<T> =
  | Readonly<{ readonly status: 'accepted'; readonly value: T; readonly commandId: ThreadCommandId }>
  | ThreadFailure;

export type CreateThreadCommand = Readonly<{
  readonly commandId: ThreadCommandId;
  readonly threadId: ThreadId;
  readonly ownerId: ThreadOwnerId;
  readonly limits: ThreadLimits;
}>;

export type AppendSettledThreadEntryCommand = Readonly<{
  readonly commandId: ThreadCommandId;
  readonly threadId: ThreadId;
  readonly entryId: ThreadEntryId;
  readonly descriptor: ThreadPayloadEnvelope;
  readonly payload: ThreadPayloadEnvelope;
}>;

export type BeginStreamingThreadEntryCommand = Readonly<{
  readonly commandId: ThreadCommandId;
  readonly threadId: ThreadId;
  readonly entryId: ThreadEntryId;
  readonly descriptor: ThreadPayloadEnvelope;
}>;

export type AppendThreadDeltaCommand = Readonly<{
  readonly commandId: ThreadCommandId;
  readonly threadId: ThreadId;
  readonly entryId: ThreadEntryId;
  readonly start: number;
  readonly payload: ThreadPayloadEnvelope;
}>;

export type HeartbeatThreadEntryCommand = Readonly<{
  readonly threadId: ThreadId;
  readonly entryId: ThreadEntryId;
}>;

export type SettleThreadEntryCommand = Readonly<{
  readonly commandId: ThreadCommandId;
  readonly threadId: ThreadId;
  readonly entryId: ThreadEntryId;
  readonly payload: ThreadPayloadEnvelope;
}>;

export type AbortThreadEntryCommand = Readonly<{
  readonly commandId: ThreadCommandId;
  readonly threadId: ThreadId;
  readonly entryId: ThreadEntryId;
  readonly reason: ThreadAbortReason;
}>;

export type ReadThreadWindowRequest = Readonly<{
  readonly threadId: ThreadId;
  readonly start?: number;
  readonly limit: number;
  readonly direction: 'forward' | 'backward';
}>;

export type ReadThreadDeltasRequest = Readonly<{
  readonly threadId: ThreadId;
  readonly entryId: ThreadEntryId;
  readonly afterOffset: number;
  readonly limit: number;
}>;

export type RemoveThreadCommand = Readonly<{
  readonly commandId: ThreadCommandId;
  readonly threadId: ThreadId;
}>;

export type ThreadRemoveResult = Readonly<{
  readonly threadId: ThreadId;
  readonly lifecycle: 'closing' | 'removed';
}>;

export type ThreadPort = Readonly<{
  create(input: CreateThreadCommand): Promise<ThreadCommandResult<ThreadSnapshot>>;
  append(input: AppendSettledThreadEntryCommand): Promise<ThreadCommandResult<ThreadWriteResult>>;
  begin(input: BeginStreamingThreadEntryCommand): Promise<ThreadCommandResult<ThreadWriteResult>>;
  appendDelta(input: AppendThreadDeltaCommand): Promise<ThreadCommandResult<ThreadAppendDeltaResult>>;
  heartbeat(input: HeartbeatThreadEntryCommand): Promise<ThreadDeltaResult | ThreadFailure>;
  settle(input: SettleThreadEntryCommand): Promise<ThreadCommandResult<ThreadWriteResult>>;
  abort(input: AbortThreadEntryCommand): Promise<ThreadCommandResult<ThreadWriteResult>>;
  /** The page, or why not: `not-found` (`thread-missing`), `invalid` (a limit out of range, an entry over the bytes). */
  window(input: ReadThreadWindowRequest): Promise<ThreadWindow | ThreadFailure>;
  /** The page, or why not: `not-found` (`thread-missing`, `entry-missing`), `invalid` (a limit out of range). */
  deltas(input: ReadThreadDeltasRequest): Promise<ThreadDeltaPage | ThreadFailure>;
  remove(input: RemoveThreadCommand): Promise<ThreadCommandResult<ThreadRemoveResult>>;
}>;
