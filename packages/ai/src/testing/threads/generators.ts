/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  createThreadLimits,
  THREAD_CONFORMANCE_LIMITS,
  THREAD_PAYLOAD_VERSION,
  type ThreadLimits,
  type ThreadPayloadEnvelope,
} from '../../threads/index.js';

export function threadConformanceEnvelope(value: unknown, codec = 'conformance.thread.v1'): ThreadPayloadEnvelope {
  return { codec, version: THREAD_PAYLOAD_VERSION, value };
}

export function threadConformanceLimits(overrides: Partial<ThreadLimits> = {}): ThreadLimits {
  return createThreadLimits({ ...THREAD_CONFORMANCE_LIMITS, ...overrides });
}
