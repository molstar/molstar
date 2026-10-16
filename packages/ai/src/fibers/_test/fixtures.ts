/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type FiberClock,
  type FiberEntryId,
  type FiberId,
  type FiberLimits,
  type FiberOwnerId,
  type FiberPayloadEnvelope,
  type FiberPort,
  fiberCommandId,
  fiberEntryId,
  fiberId,
  fiberOwnerId,
} from '../index.js';
import { fiberConformanceEnvelope, fiberConformanceLimits } from './generators.js';

export class FiberManualClock implements FiberClock {
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

/** A fiber port over fresh, empty storage; `reopen` gives another over the same storage (a restart). */
export type FiberBackend = {
  readonly fibers: FiberPort;
  /** A new port over the storage the first one wrote, as after a restart. Absent when the storage does not outlive its port. */
  reopen?(): FiberPort;
};

export type DeterministicFiberFixture = {
  readonly fibers: FiberPort;
  /** A port over the same storage as `fibers`, as after a restart; absent when the storage does not outlive it. */
  readonly reopen?: () => FiberPort;
  readonly clock: FiberManualClock;
  readonly ownerId: FiberOwnerId;
  readonly fiberId: FiberId;
  readonly otherFiberId: FiberId;
  readonly limits: FiberLimits;
  envelope(value: unknown): FiberPayloadEnvelope;
  commandId(value: string): ReturnType<typeof fiberCommandId>;
  entryId(value: string): FiberEntryId;
};

export type FiberFixtureFactoryOptions = {
  readonly limits?: FiberLimits;
  readonly clock?: FiberManualClock;
};

export type FiberFixtureDependencies = {
  readonly createFiber: (options: { readonly clock: FiberClock; readonly limits: FiberLimits }) => FiberBackend;
};

export function createDeterministicFiberFixture(
  dependencies: FiberFixtureDependencies,
  options: FiberFixtureFactoryOptions = {},
): DeterministicFiberFixture {
  const clock = options.clock ?? new FiberManualClock();
  const limits = options.limits ?? fiberConformanceLimits();
  const backend = dependencies.createFiber({ clock, limits });
  const reopen = backend.reopen?.bind(backend);
  return {
    fibers: backend.fibers,
    ...(reopen === undefined ? {} : { reopen }),
    clock,
    ownerId: fiberOwnerId('fiber-conformance-owner'),
    fiberId: fiberId('fiber-conformance-fiber'),
    otherFiberId: fiberId('fiber-conformance-other-fiber'),
    limits,
    envelope: fiberConformanceEnvelope,
    commandId: fiberCommandId,
    entryId: fiberEntryId,
  };
}
