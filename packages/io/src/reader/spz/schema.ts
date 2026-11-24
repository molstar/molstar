/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

/**
 * Niantic's gaussian splats (`.spz`), versions 1 to 3: a 16-byte header and arrays of all splats, each a view of the
 * (inflated) bytes, as stored (no copy, not decoded):
 * - `positions`: x, y, z a splat, 24-bit signed fixed point of `fractionalBits` (little-endian), version 1 `float16`
 * - `alphas`: the opacity (sigmoid's) as a byte
 * - `colors`: r, g, b, the zeroth spherical harmonic as bytes `dc * 0.15 * 255 + 127.5`
 * - `scales`: x, y, z, the natural logs as bytes `(log + 10) * 16`
 * - `rotations`: version 2 x, y, z as bytes `q * 127.5 + 127.5` (w the rest, not negative); version 3 a little-endian
 *   `u32`, the smallest three: the largest component's index (0 x to 3 w) in bits 30 and 31, the others, in their
 *   order, 10 bits each from the high to the low (a sign in bit 9 and the magnitude over `sqrt(1/2)` of 9 bits)
 * - `sh`: the spherical harmonics past the zeroth, `SpzShDims[shDegree]` coefficients a channel, coefficient-major
 *   (r, g, b of each), as bytes `x * 128 + 128`
 *
 * Stored right-up-back (the PLY's right-down-forward turned 180° about x)
 */
export interface SpzFile {
  readonly version: number;
  readonly count: number;
  readonly shDegree: number;
  readonly fractionalBits: number;
  /** 1: trained with antialiasing */
  readonly flags: number;
  readonly positions: ReadonlyUint8Array;
  readonly alphas: ReadonlyUint8Array;
  readonly colors: ReadonlyUint8Array;
  readonly scales: ReadonlyUint8Array;
  readonly rotations: ReadonlyUint8Array;
  readonly sh: ReadonlyUint8Array;
}

/** "NGSP", little-endian */
export const SpzMagic = 0x5053474e;

/** Coefficients past the zeroth a channel, by degree (0 to 4) */
export const SpzShDims = [0, 3, 8, 15, 24] as const;

/** The scale of a color byte to the zeroth spherical harmonic (not its constant 0.282) */
export const SpzColorScale = 0.15;
