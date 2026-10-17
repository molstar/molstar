/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Seq*'s guard harness for Prism's tiers (packages/prism/AGENTS.md, Tests): `quickCase` for the quick allocation tier
// (bytes a step against `perf-baseline.json`, its positive control built in) and the heap guards' scaling (n and 4n
// steps, `leaking` for the control). Test code only: the harness is Seq*'s `_test/model` (packages/sequences/AGENTS.md).

export { quickCase } from '../../../../sequences/core/_test/model/quick.js';
export {
  controlCaught,
  describeScaling,
  heapScaling,
  isFlat,
  leaking,
} from '../../../../sequences/core/_test/model/heap.js';
