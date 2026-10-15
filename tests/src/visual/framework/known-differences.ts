/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type Tolerance } from './compare.js';
import { matchesAny } from './test-case.js';

export type KnownDifferenceAction =
  | { kind: 'tolerance'; tolerance: Partial<Tolerance> }
  /** Rendered and shown, but a mismatch does not fail */
  | { kind: 'expected-diff' }
  | { kind: 'unsupported' };

export interface KnownDifference {
  /** Target profile id globs, e.g. `webgl1*` */
  readonly profiles: readonly string[];
  /** Test id globs */
  readonly tests: readonly string[];
  readonly reason: string;
  /** Anchor in the backend-differences documentation */
  readonly doc?: string;
  readonly action: KnownDifferenceAction;
}

/**
 * Source of truth for documented rendering differences between backends,
 * relative to the reference profile.
 */
export const KnownDifferences: KnownDifference[] = [
  {
    profiles: ['webgl1*'],
    tests: ['renderer/illumination/**', 'fallback/illumination'],
    reason: 'WebGL1 trace shader uses a float hash instead of the uint-based RNG, so the noise pattern differs.',
    doc: 'illumination-rng',
    action: { kind: 'tolerance', tolerance: { maxDiffRatio: 0.03 } },
  },
  {
    profiles: ['webgpu'],
    tests: ['postprocessing/occlusion/multi-scale-half-resolution'],
    reason:
      'Half-resolution SSAO is upsampled with bilinear filtering of occlusion packed into two 8-bit channels; with ' +
      'the flipped row order the sampler weights are mirrored and quantize differently, which the packing amplifies.',
    action: { kind: 'tolerance', tolerance: { maxDiffRatio: 0.005 } },
  },
  {
    profiles: ['webgpu'],
    tests: ['geometry/cylinders/star'],
    reason:
      'All eight cylinders overlap at the shared origin, where the front surface is decided by nearly equal fragment ' +
      'depths; the [0, 1] depth projection rounds them differently than the GL [-1, 1] one (subsets of the star match).',
    action: { kind: 'tolerance', tolerance: { maxDiffRatio: 0.006 } },
  },
];

export function findKnownDifference(profileId: string, testId: string): KnownDifference | undefined {
  return KnownDifferences.find((d) => matchesAny(profileId, d.profiles) && matchesAny(testId, d.tests));
}
