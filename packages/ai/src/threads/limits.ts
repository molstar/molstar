/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type ThreadLimits = Readonly<{
  readonly maxEntryBytes: number;
  readonly maxBatchBytes: number;
  readonly deltaFlushBytes: number;
  readonly deltaFlushDelayMs: number;
  readonly maxActiveStreams: number;
  readonly maxActiveStreamBytes: number;
  readonly maxActiveStreamBatches: number;
  readonly maxWindowEntries: number;
  readonly maxWindowBytes: number;
  readonly maxDeltaPageBatches: number;
  readonly maxDeltaPageBytes: number;
  readonly heartbeatIntervalMs: number;
  readonly inactivityTimeoutMs: number;
  readonly abortedPreviewRetentionMs: number;
}>;

/** Default ceiling used when a runtime does not pass its own bounds. */
export const THREAD_LIMIT_HARD_BOUNDS: ThreadLimits = Object.freeze({
  maxEntryBytes: 1_000_000,
  maxBatchBytes: 256_000,
  deltaFlushBytes: 256_000,
  deltaFlushDelayMs: 60_000,
  maxActiveStreams: 32,
  maxActiveStreamBytes: 2_000_000,
  maxActiveStreamBatches: 10_000,
  maxWindowEntries: 256,
  maxWindowBytes: 2_000_000,
  maxDeltaPageBatches: 256,
  maxDeltaPageBytes: 512_000,
  heartbeatIntervalMs: 60_000,
  inactivityTimeoutMs: 3_600_000,
  abortedPreviewRetentionMs: 3_600_000,
});

/** Small portable profile for conformance. Not a production default. */
export const THREAD_CONFORMANCE_LIMITS: ThreadLimits = Object.freeze({
  maxEntryBytes: 4_096,
  maxBatchBytes: 512,
  deltaFlushBytes: 256,
  deltaFlushDelayMs: 20,
  maxActiveStreams: 4,
  maxActiveStreamBytes: 2_048,
  maxActiveStreamBatches: 8,
  maxWindowEntries: 8,
  maxWindowBytes: 8_192,
  maxDeltaPageBatches: 4,
  maxDeltaPageBytes: 1_024,
  heartbeatIntervalMs: 50,
  inactivityTimeoutMs: 200,
  abortedPreviewRetentionMs: 0,
});

const POSITIVE_LIMIT_KEYS = [
  'maxEntryBytes',
  'maxBatchBytes',
  'deltaFlushBytes',
  'deltaFlushDelayMs',
  'maxActiveStreams',
  'maxActiveStreamBytes',
  'maxActiveStreamBatches',
  'maxWindowEntries',
  'maxWindowBytes',
  'maxDeltaPageBatches',
  'maxDeltaPageBytes',
  'heartbeatIntervalMs',
  'inactivityTimeoutMs',
] as const satisfies readonly (keyof ThreadLimits)[];

function assertIntegerInRange(label: string, value: number, min: number, max: number): void {
  if (!Number.isInteger(value) || value < min || value > max) {
    throw new Error(`Thread limit ${label} must be an integer in ${min}..${max}`);
  }
}

function assertThreadLimitRelations(limits: ThreadLimits): void {
  if (limits.heartbeatIntervalMs >= limits.inactivityTimeoutMs) {
    throw new Error('Thread heartbeatIntervalMs must be strictly less than inactivityTimeoutMs');
  }
  if (limits.maxBatchBytes > limits.maxEntryBytes) {
    throw new Error('Thread maxBatchBytes must not exceed maxEntryBytes');
  }
  if (limits.maxWindowBytes < 2 * limits.maxEntryBytes) {
    throw new Error('Thread maxWindowBytes must be at least twice maxEntryBytes');
  }
  if (limits.maxDeltaPageBytes < limits.maxBatchBytes) {
    throw new Error('Thread maxDeltaPageBytes must be at least maxBatchBytes');
  }
  if (limits.maxActiveStreamBytes < limits.maxBatchBytes) {
    throw new Error('Thread maxActiveStreamBytes must be at least maxBatchBytes');
  }
}

function assertThreadLimitShape(limits: ThreadLimits, ceiling: ThreadLimits): void {
  for (const key of POSITIVE_LIMIT_KEYS) {
    assertIntegerInRange(key, limits[key], 1, ceiling[key]);
  }
  assertIntegerInRange(
    'abortedPreviewRetentionMs',
    limits.abortedPreviewRetentionMs,
    0,
    ceiling.abortedPreviewRetentionMs,
  );
  assertThreadLimitRelations(limits);
}

export function createThreadLimits(
  limits: ThreadLimits,
  bounds: ThreadLimits = THREAD_LIMIT_HARD_BOUNDS,
): ThreadLimits {
  assertThreadLimitShape(bounds, bounds);
  assertThreadLimitShape(limits, bounds);
  return Object.freeze({ ...limits });
}
