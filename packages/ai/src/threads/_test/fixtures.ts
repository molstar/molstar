/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Clock } from '../../messages/index.js';
import {
  type ThreadEntryId,
  type ThreadId,
  type ThreadLimits,
  type ThreadOwnerId,
  type ThreadPayloadEnvelope,
  type ThreadPort,
  threadCommandId,
  threadEntryId,
  threadId,
  threadOwnerId,
} from '../index.js';
import { threadConformanceEnvelope, threadConformanceLimits } from './generators.js';

export class ThreadManualClock implements Pick<Clock, 'now'> {
  #nowMs: number;

  constructor(startMs = 1_000) {
    this.#nowMs = startMs;
  }

  now(): number {
    return this.#nowMs;
  }

  set(value: number): void {
    this.#nowMs = value;
  }

  advance(deltaMs: number): void {
    if (!Number.isSafeInteger(deltaMs) || deltaMs < 0) throw new Error('clock advance must be nonnegative');
    this.#nowMs += deltaMs;
  }
}

/** A thread port over fresh, empty storage; `reopen` gives another over the same storage (a restart). */
export type ThreadBackend = {
  readonly threads: ThreadPort;
  /** A new port over the storage the first one wrote, as after a restart. Absent when the storage does not outlive its port. */
  reopen?(): ThreadPort;
};

export type DeterministicThreadFixture = {
  readonly threads: ThreadPort;
  /** A port over the same storage as `threads`, as after a restart; absent when the storage does not outlive it. */
  readonly reopen?: () => ThreadPort;
  readonly clock: ThreadManualClock;
  readonly ownerId: ThreadOwnerId;
  readonly threadId: ThreadId;
  readonly otherThreadId: ThreadId;
  readonly limits: ThreadLimits;
  envelope(value: unknown): ThreadPayloadEnvelope;
  commandId(value: string): ReturnType<typeof threadCommandId>;
  entryId(value: string): ThreadEntryId;
};

export type ThreadFixtureFactoryOptions = {
  readonly limits?: ThreadLimits;
  readonly clock?: ThreadManualClock;
};

export type ThreadFixtureDependencies = {
  readonly createThread: (options: {
    readonly clock: Pick<Clock, 'now'>;
    readonly limits: ThreadLimits;
  }) => ThreadBackend;
};

export function createDeterministicThreadFixture(
  dependencies: ThreadFixtureDependencies,
  options: ThreadFixtureFactoryOptions = {},
): DeterministicThreadFixture {
  const clock = options.clock ?? new ThreadManualClock();
  const limits = options.limits ?? threadConformanceLimits();
  const backend = dependencies.createThread({ clock, limits });
  const reopen = backend.reopen?.bind(backend);
  return {
    threads: backend.threads,
    ...(reopen === undefined ? {} : { reopen }),
    clock,
    ownerId: threadOwnerId('thread-conformance-owner'),
    threadId: threadId('thread-conformance-thread'),
    otherThreadId: threadId('thread-conformance-other-thread'),
    limits,
    envelope: threadConformanceEnvelope,
    commandId: threadCommandId,
    entryId: threadEntryId,
  };
}
