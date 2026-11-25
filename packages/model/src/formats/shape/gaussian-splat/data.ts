/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Gaussian splats of a file (PLY, `.splat`, `.spz`, `.sog`) as the graphics build them (`GaussianSplats.fromData`): each
 * splat decoded on demand from the file's own data (no copy of the file's arrays), into a `GaussianSplat`, in the
 * file's frame (the graphics turn it upright, see `upright`)
 */
export interface GaussianSplatData {
  readonly count: number;
  /**
   * Stored right-up-back (SPZ), else right-down-forward (the 3DGS PLY's y-down, the camera convention of training):
   * turned 180° about x by the shape's transform where it is not as wanted
   */
  readonly upright: boolean;
  /** splat `i` into `out` */
  getSplat(i: number, out: GaussianSplat): void;
  /** the view-dependent color, none without */
  readonly sh?: GaussianSplatSh;
}

/** A splat: center, natural logs of its scales, rotation (a unit quaternion `x, y, z, w`), opacity, base color */
export interface GaussianSplat {
  x: number;
  y: number;
  z: number;
  logScaleX: number;
  logScaleY: number;
  logScaleZ: number;
  qx: number;
  qy: number;
  qz: number;
  qw: number;
  /** 0 to 1 */
  opacity: number;
  /** the base color, 0 to 1 (not clamped) */
  r: number;
  g: number;
  b: number;
}

export function GaussianSplat(): GaussianSplat {
  return {
    x: 0,
    y: 0,
    z: 0,
    logScaleX: 0,
    logScaleY: 0,
    logScaleZ: 0,
    qx: 0,
    qy: 0,
    qz: 0,
    qw: 1,
    opacity: 1,
    r: 1,
    g: 1,
    b: 1,
  };
}

/** View-dependent color: spherical harmonics past the zeroth, the basis of 3D Gaussian splatting (its `eval_sh`) */
export interface GaussianSplatSh {
  readonly degree: 1 | 2 | 3;
  /** coefficient of basis function `basis` (past the zeroth, below `ShBasisCounts[degree]`) of channel `channel` */
  get(splat: number, basis: number, channel: number): number;
}

/** Options of reading gaussian splats */
export interface GaussianSplatReadOptions {
  /** the degree of view-dependent color read, at most the file's (0: none) */
  readonly shDegree: number;
}

/** The zeroth spherical harmonic's constant (Y_0^0): a splat's base color is `0.5 + SH_C0 * f_dc` */
export const SH_C0 = 0.28209479177387814;

/** Spherical harmonics basis functions a channel past the zeroth, by degree 0 to 3 */
export const ShBasisCounts = [0, 3, 8, 15] as const;

export function sigmoid(x: number) {
  return 1 / (1 + Math.exp(-x));
}

/** The degree of `count` coefficients past the zeroth (3 channels): 1, 2 or 3, else 0 */
export function getShDegree(count: number): 0 | 1 | 2 | 3 {
  const n = count / 3;
  return n === 15 ? 3 : n === 8 ? 2 : n === 3 ? 1 : 0;
}

/** Sets the center and rotation `(x, y, z, w)` of `out` (the rotation as given: normalized where read) */
export function setCenterRotation(
  out: GaussianSplat,
  x: number,
  y: number,
  z: number,
  qx: number,
  qy: number,
  qz: number,
  qw: number,
) {
  out.x = x;
  out.y = y;
  out.z = z;
  out.qx = qx;
  out.qy = qy;
  out.qz = qz;
  out.qw = qw;
}

/**
 * The view-dependent color up to `degree` (at most 3, below 1: none) read by `get` (the file's coefficient of a splat,
 * basis function and channel)
 */
export function createShAccessor(
  degree: number,
  get: (splat: number, basis: number, channel: number) => number,
): GaussianSplatSh | undefined {
  if (degree < 1) return;
  return { degree: Math.min(3, degree) as 1 | 2 | 3, get };
}
