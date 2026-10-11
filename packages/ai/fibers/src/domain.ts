/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { FiberCommandId, FiberEntryId, FiberId, FiberOwnerId } from './ids.js';
import type { FiberLimits } from './limits.js';

export const FIBER_PAYLOAD_VERSION = 1;

export type FiberPayloadEnvelope = Readonly<{
  readonly codec: string;
  readonly version: typeof FIBER_PAYLOAD_VERSION;
  readonly value: unknown;
}>;

export type FiberAbortReason = 'caller-interrupted' | 'inactivity' | 'limit-exceeded' | 'host-closed';

export type FiberEntryState =
  | Readonly<{
      readonly kind: 'streaming';
      readonly nextOffset: number;
      readonly lastActivityAt: number;
    }>
  | Readonly<{
      readonly kind: 'settled';
      readonly payload: FiberPayloadEnvelope;
      readonly settledAt: number;
    }>
  | Readonly<{
      readonly kind: 'aborted';
      readonly reason: FiberAbortReason;
      readonly endedAt: number;
    }>;

export type FiberEntry = Readonly<{
  readonly entryId: FiberEntryId;
  readonly fiberId: FiberId;
  readonly position: number;
  readonly descriptor: FiberPayloadEnvelope;
  readonly state: FiberEntryState;
}>;

export type FiberDeltaBatch = Readonly<{
  readonly entryId: FiberEntryId;
  readonly start: number;
  readonly endExclusive: number;
  readonly payload: FiberPayloadEnvelope;
}>;

export type FiberSnapshot = Readonly<{
  readonly fiberId: FiberId;
  readonly ownerId: FiberOwnerId;
  readonly nextPosition: number;
  readonly lifecycle: 'active' | 'closing' | 'removed';
}>;

export type FiberWriteResult = Readonly<{
  readonly fiber: FiberSnapshot;
  readonly entry: FiberEntry;
}>;

export type FiberDeltaResult = Readonly<{
  readonly fiber: FiberSnapshot;
  readonly entry: FiberEntry;
  readonly nextOffset: number;
}>;

export type FiberAppendDeltaResult =
  | Readonly<{
      readonly kind: 'appended';
      readonly fiber: FiberSnapshot;
      readonly entry: FiberEntry;
      readonly nextOffset: number;
    }>
  | Readonly<{
      readonly kind: 'aborted';
      readonly fiber: FiberSnapshot;
      readonly entry: FiberEntry;
      readonly reason: Extract<FiberAbortReason, 'limit-exceeded'>;
    }>;

export type FiberWindow = Readonly<{
  readonly fiberId: FiberId;
  readonly entries: readonly FiberEntry[];
  readonly start: number;
  readonly endExclusive: number;
  readonly hasMore: boolean;
}>;

export type FiberDeltaPage = Readonly<{
  readonly entryId: FiberEntryId;
  readonly batches: readonly FiberDeltaBatch[];
  readonly nextOffset: number;
  readonly hasMore: boolean;
  readonly gap?: Readonly<{ readonly afterOffset: number }>;
}>;

export type FiberCleanupKind = 'deltas' | 'entries' | 'streams' | 'receipts' | 'fiber';

export type FiberCleanupCursor = Readonly<{
  readonly kind: FiberCleanupKind;
  readonly after?: string;
}>;

export type FiberCleanupPage = Readonly<{
  readonly kind: FiberCleanupKind;
  readonly ids: readonly string[];
  readonly next?: FiberCleanupCursor;
}>;

export type FiberFailure = Readonly<{
  readonly status: 'conflict' | 'invalid' | 'not-found' | 'closed';
  readonly reason: string;
}>;

export type FiberCommandResult<T> =
  | Readonly<{ readonly status: 'accepted'; readonly value: T; readonly commandId: FiberCommandId }>
  | FiberFailure;

export type FiberClock = Readonly<{ now(): number }>;

export type CreateFiberCommand = Readonly<{
  readonly commandId: FiberCommandId;
  readonly fiberId: FiberId;
  readonly ownerId: FiberOwnerId;
  readonly limits: FiberLimits;
}>;

export type AppendSettledFiberEntryCommand = Readonly<{
  readonly commandId: FiberCommandId;
  readonly fiberId: FiberId;
  readonly entryId: FiberEntryId;
  readonly descriptor: FiberPayloadEnvelope;
  readonly payload: FiberPayloadEnvelope;
}>;

export type BeginStreamingFiberEntryCommand = Readonly<{
  readonly commandId: FiberCommandId;
  readonly fiberId: FiberId;
  readonly entryId: FiberEntryId;
  readonly descriptor: FiberPayloadEnvelope;
}>;

export type AppendFiberDeltaCommand = Readonly<{
  readonly commandId: FiberCommandId;
  readonly fiberId: FiberId;
  readonly entryId: FiberEntryId;
  readonly start: number;
  readonly payload: FiberPayloadEnvelope;
}>;

export type HeartbeatFiberEntryCommand = Readonly<{
  readonly fiberId: FiberId;
  readonly entryId: FiberEntryId;
}>;

export type SettleFiberEntryCommand = Readonly<{
  readonly commandId: FiberCommandId;
  readonly fiberId: FiberId;
  readonly entryId: FiberEntryId;
  readonly payload: FiberPayloadEnvelope;
}>;

export type AbortFiberEntryCommand = Readonly<{
  readonly commandId: FiberCommandId;
  readonly fiberId: FiberId;
  readonly entryId: FiberEntryId;
  readonly reason: FiberAbortReason;
}>;

export type ReadFiberWindowRequest = Readonly<{
  readonly fiberId: FiberId;
  readonly start?: number;
  readonly limit: number;
  readonly direction: 'forward' | 'backward';
}>;

export type ReadFiberDeltasRequest = Readonly<{
  readonly fiberId: FiberId;
  readonly entryId: FiberEntryId;
  readonly afterOffset: number;
  readonly limit: number;
}>;

export type RemoveFiberCommand = Readonly<{
  readonly commandId: FiberCommandId;
  readonly fiberId: FiberId;
}>;

export type FiberRemoveResult = Readonly<{
  readonly fiberId: FiberId;
  readonly lifecycle: 'closing' | 'removed';
}>;

export type FiberPort = Readonly<{
  create(input: CreateFiberCommand): Promise<FiberCommandResult<FiberSnapshot>>;
  append(input: AppendSettledFiberEntryCommand): Promise<FiberCommandResult<FiberWriteResult>>;
  begin(input: BeginStreamingFiberEntryCommand): Promise<FiberCommandResult<FiberWriteResult>>;
  appendDelta(input: AppendFiberDeltaCommand): Promise<FiberCommandResult<FiberAppendDeltaResult>>;
  heartbeat(input: HeartbeatFiberEntryCommand): Promise<FiberDeltaResult | FiberFailure>;
  settle(input: SettleFiberEntryCommand): Promise<FiberCommandResult<FiberWriteResult>>;
  abort(input: AbortFiberEntryCommand): Promise<FiberCommandResult<FiberWriteResult>>;
  window(input: ReadFiberWindowRequest): Promise<FiberWindow | undefined>;
  deltas(input: ReadFiberDeltasRequest): Promise<FiberDeltaPage | undefined>;
  remove(input: RemoveFiberCommand): Promise<FiberCommandResult<FiberRemoveResult>>;
}>;
