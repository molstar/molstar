/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export { cloneData } from './data.js';
export type {
  AbortThreadEntryCommand,
  AppendThreadDeltaCommand,
  AppendSettledThreadEntryCommand,
  BeginStreamingThreadEntryCommand,
  CreateThreadCommand,
  ThreadAbortReason,
  ThreadAppendDeltaResult,
  ThreadClock,
  ThreadCommandResult,
  ThreadDeltaBatch,
  ThreadDeltaPage,
  ThreadDeltaResult,
  ThreadEntry,
  ThreadEntryState,
  ThreadFailure,
  ThreadPayloadEnvelope,
  ThreadPort,
  ThreadRemoveResult,
  ThreadSnapshot,
  ThreadWindow,
  ThreadWriteResult,
  HeartbeatThreadEntryCommand,
  ReadThreadDeltasRequest,
  ReadThreadWindowRequest,
  RemoveThreadCommand,
  SettleThreadEntryCommand,
} from './domain.js';
export { THREAD_PAYLOAD_VERSION } from './domain.js';
export { createThreadHost, type ThreadHost, type ThreadHostOptions } from './host.js';
export {
  type ThreadCommandId,
  type ThreadEntryId,
  type ThreadId,
  type ThreadOwnerId,
  threadCommandId,
  threadEntryId,
  threadId,
  threadOwnerId,
} from './ids.js';
export {
  type ThreadCommandKind,
  type ThreadKernelDecision,
  type ThreadKernelThread,
  type ThreadKernelReceipt,
  type ThreadKernelStream,
  type ThreadKernelTouchDecision,
  threadCommandFingerprint,
  threadPageTakeBound,
  threadStoredEnvelopeBytes,
  threadWindowRangeStart,
  prepareThreadAbort,
  prepareThreadAppendDelta,
  prepareThreadAppendSettled,
  prepareThreadBegin,
  prepareThreadCreate,
  prepareThreadHeartbeat,
  prepareThreadInactivityAbort,
  prepareThreadRemove,
  prepareThreadSettle,
  projectThreadDeltas,
  projectThreadWindow,
  validThreadPayloadEnvelope,
} from './kernel.js';
export {
  createThreadLimits,
  THREAD_CONFORMANCE_LIMITS,
  THREAD_LIMIT_HARD_BOUNDS,
  type ThreadLimits,
} from './limits.js';
export { createMemoryRows } from './memory.js';
export {
  type Field,
  type ThreadRows,
  type ThreadTable,
  type Row,
  RowConflict,
  type RowStore,
  type RowWrite,
  type Unsubscribe,
} from './storage.js';
