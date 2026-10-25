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
    profiles: ['webgpu*'],
    tests: ['geometry/palette/*'],
    reason:
      'The palette coordinate is interpolated across triangles and looked up with nearest filtering: where it falls on ' +
      'a texel boundary, rasterization rounding picks the neighboring palette color, a one-pixel line along the band edges.',
    action: { kind: 'tolerance', tolerance: { maxDiffRatio: 0.008 } },
  },
  {
    profiles: ['*'],
    tests: ['transparency/*/solid-interior'],
    reason:
      'Not reproducible on GPUs: with ANGLE on D3D11, pages rendering in parallel sometimes compute slightly different ' +
      'cap depths, which the exact depth comparisons of the peeling turn into a few hundred differing pixels; ' +
      'sequential runs and software rendering are deterministic (see scripts/visual-tests.mjs).',
    action: { kind: 'tolerance', tolerance: { maxDiffRatio: 0.005 } },
  },
  {
    profiles: ['webgl1*'],
    tests: ['renderer/illumination/**', 'fallback/illumination'],
    reason: 'WebGL1 trace shader uses a float hash instead of the uint-based RNG, so the noise pattern differs.',
    doc: 'illumination-rng',
    action: { kind: 'tolerance', tolerance: { maxDiffRatio: 0.03 } },
  },
  {
    profiles: ['webgpu*'],
    tests: ['geometry/direct-volume/opaque-transfer'],
    reason:
      'A transfer function stepping from 0 to 1 within one cell: where a sample falls on the step, the backends round ' +
      'the interpolated value to either side, so the ray turns opaque one step earlier or later.',
    action: { kind: 'tolerance', tolerance: { maxDiffRatio: 0.004 } },
  },
];

export function findKnownDifference(profileId: string, testId: string): KnownDifference | undefined {
  return KnownDifferences.find((d) => matchesAny(profileId, d.profiles) && matchesAny(testId, d.tests));
}
