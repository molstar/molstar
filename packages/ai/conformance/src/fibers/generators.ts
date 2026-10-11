/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  createFiberLimits,
  FIBER_CONFORMANCE_LIMITS,
  FIBER_PAYLOAD_VERSION,
  type FiberCleanupCursor,
  type FiberLimits,
  type FiberPayloadEnvelope,
  validFiberPayloadEnvelope,
} from '@molstar/ai-fibers';

export function fiberConformanceEnvelope(value: unknown, codec = 'conformance.fiber.v1'): FiberPayloadEnvelope {
  return { codec, version: FIBER_PAYLOAD_VERSION, value };
}

export const FIBER_INVALID_ENVELOPES: readonly unknown[] = Object.freeze([
  { codec: 'x', version: 2, value: null },
  { codec: '', version: 1, value: null },
  { codec: 'x', version: 1, value: undefined },
  { codec: 'x', version: 1, value: null, extra: true },
  'not-an-envelope',
]);

export function fiberConformanceLimits(overrides: Partial<FiberLimits> = {}): FiberLimits {
  return createFiberLimits({ ...FIBER_CONFORMANCE_LIMITS, ...overrides });
}

export function fiberCleanupCursorPages(ids: readonly string[], pageSize: number): readonly FiberCleanupCursor[] {
  const cursors: FiberCleanupCursor[] = [{ kind: 'entries' }];
  const ordered = [...ids].sort();
  for (let index = pageSize; index < ordered.length; index += pageSize) {
    const after = ordered[index - 1];
    if (after !== undefined) cursors.push({ kind: 'entries', after });
  }
  return cursors;
}

export function fiberEnvelopeIsValid(value: unknown, maxBytes = FIBER_CONFORMANCE_LIMITS.maxEntryBytes): boolean {
  return validFiberPayloadEnvelope(value, maxBytes);
}
