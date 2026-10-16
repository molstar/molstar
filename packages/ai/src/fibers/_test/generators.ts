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
  type FiberLimits,
  type FiberPayloadEnvelope,
} from '../index.js';

export function fiberConformanceEnvelope(value: unknown, codec = 'conformance.fiber.v1'): FiberPayloadEnvelope {
  return { codec, version: FIBER_PAYLOAD_VERSION, value };
}

export function fiberConformanceLimits(overrides: Partial<FiberLimits> = {}): FiberLimits {
  return createFiberLimits({ ...FIBER_CONFORMANCE_LIMITS, ...overrides });
}
