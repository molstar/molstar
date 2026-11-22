/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat3, Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { ShBasisCounts } from '@molstar/model/formats/shape/gaussian-splat/data';

/**
 * A splat's shape and opacity as the shaders read them (`tSplat`, `rgba16i`): two elements of four int16, the rotation
 * (x, y, z snorm16 and w unorm16: w kept non-negative, `q` and `-q` the same rotation, so the unsigned fourth holds it)
 * and the scales (x, y, z the natural logs in fixed point, `SplatLogScaleFactor` steps a unit) with the opacity
 * (unorm16). 16 bytes a splat; the covariance `R S Sᵀ Rᵀ` built where it is read
 */
export const SplatTexels = 2;
/** Int16 per splat: four per element */
export const SplatStride = SplatTexels * 4;

/**
 * Fixed-point steps per unit of a natural-log scale: logs within ±32 (scales from 1e-14 to 8e13 of the positions'
 * unit) in relative steps of 0.1%
 */
export const SplatLogScaleFactor = 1024;

const Int16Max = 32767;
const Uint16Max = 65535;

function clampInt16(v: number) {
  // NaN and -Infinity (the log of a zero scale) the smallest
  return v >= -Int16Max ? Math.min(Int16Max, v) : -Int16Max;
}

/**
 * Sets splat `i` of `data` (`SplatStride` int16 each): the natural logs of its scales (as 3DGS files store them), its
 * rotation (any length, normalized; zero the identity) and its opacity (0 to 1)
 */
export function setSplat(
  data: Int16Array,
  i: number,
  logScaleX: number,
  logScaleY: number,
  logScaleZ: number,
  qx: number,
  qy: number,
  qz: number,
  qw: number,
  opacity: number,
) {
  const length = Math.sqrt(qx * qx + qy * qy + qz * qz + qw * qw);
  if (length > 0 && Number.isFinite(length)) {
    // w non-negative: -q is the same rotation
    const s = (qw < 0 ? -1 : 1) / length;
    qx *= s;
    qy *= s;
    qz *= s;
    qw *= s;
  } else {
    qx = qy = qz = 0;
    qw = 1;
  }
  const o = i * SplatStride;
  data[o] = Math.round(qx * Int16Max);
  data[o + 1] = Math.round(qy * Int16Max);
  data[o + 2] = Math.round(qz * Int16Max);
  // a uint16 as its bits (an int16 array wraps it)
  data[o + 3] = Math.round(qw * Uint16Max);
  data[o + 4] = clampInt16(Math.round(logScaleX * SplatLogScaleFactor));
  data[o + 5] = clampInt16(Math.round(logScaleY * SplatLogScaleFactor));
  data[o + 6] = clampInt16(Math.round(logScaleZ * SplatLogScaleFactor));
  data[o + 7] = Math.round(Math.min(1, Math.max(0, opacity || 0)) * Uint16Max);
}

/** The rotation of splat `i` of `data` into `out`, unit length (as the shaders read it) */
export function getSplatRotation(out: Quat, data: Int16Array, i: number) {
  const o = i * SplatStride;
  Quat.set(out, data[o] / Int16Max, data[o + 1] / Int16Max, data[o + 2] / Int16Max, (data[o + 3] & 0xffff) / Uint16Max);
  return Quat.normalize(out, out);
}

/** The scales of splat `i` of `data` into `out` (standard deviations along its axes, in the positions' unit) */
export function getSplatScale(out: Vec3, data: Int16Array, i: number) {
  const o = i * SplatStride;
  return Vec3.set(
    out,
    Math.exp(data[o + 4] / SplatLogScaleFactor),
    Math.exp(data[o + 5] / SplatLogScaleFactor),
    Math.exp(data[o + 6] / SplatLogScaleFactor),
  );
}

/** The opacity of splat `i` of `data` (0 to 1) */
export function getSplatOpacity(data: Int16Array, i: number) {
  return (data[i * SplatStride + 7] & 0xffff) / Uint16Max;
}

/** The largest scale of the `count` splats of `data` (0 without splats): the bounds' padding */
export function getMaxSplatScale(data: Int16Array, count: number) {
  if (count === 0) return 0;
  let max = -Int16Max;
  for (let o = 4, ol = count * SplatStride; o < ol; o += SplatStride) {
    if (data[o] > max) max = data[o];
    if (data[o + 1] > max) max = data[o + 1];
    if (data[o + 2] > max) max = data[o + 2];
  }
  return Math.exp(max / SplatLogScaleFactor);
}

const tmpRotation = Quat();
const tmpScale = Vec3();
const tmpMat4 = Mat4();
const tmpM = Mat3();
const tmpMt = Mat3();

/**
 * The covariance `Σ = R S Sᵀ Rᵀ` of splat `i` of `data` into `out` (a symmetric 3x3 matrix), as the shaders build it
 * from what is stored
 */
export function getSplatCovariance(out: Mat3, data: Int16Array, i: number) {
  getSplatRotation(tmpRotation, data, i);
  getSplatScale(tmpScale, data, i);
  Mat3.fromMat4(tmpM, Mat4.fromQuat(tmpMat4, tmpRotation));
  // M = R S: the rotation's columns scaled
  for (let c = 0; c < 3; ++c) {
    for (let r = 0; r < 3; ++r) tmpM[c * 3 + r] *= tmpScale[c];
  }
  return Mat3.mul(out, tmpM, Mat3.transpose(tmpMt, tmpM));
}

/**
 * The view-dependent color of splats: per splat a texel of `rgba8` per basis function past the zeroth (rgb its
 * channels' coefficients, `ShBasisCounts[degree]` texels), each value `(byte / 255 * 2 - 1) * scale` (`scale` the
 * largest coefficient's magnitude); the basis of 3D Gaussian splatting (its `eval_sh`), added to the base color for the
 * direction from the camera to the splat
 */
export interface SplatSh {
  readonly degree: 0 | 1 | 2 | 3;
  readonly scale: number;
  readonly data: Uint8Array;
}

/** No view-dependent color (degree 0) */
export const NoSplatSh: SplatSh = { degree: 0, scale: 0, data: new Uint8Array(0) };

/**
 * Encodes `coefficients` (per splat `count * 3` values: per basis function its r, g, b, as `getCoefficient(i, k, c)`
 * gives them) of `degree` for `splatCount` splats
 */
export function encodeSplatSh(
  splatCount: number,
  degree: 0 | 1 | 2 | 3,
  getCoefficient: (splat: number, basis: number, channel: number) => number,
): SplatSh {
  const n = ShBasisCounts[degree];
  let scale = 0;
  for (let i = 0; i < splatCount; ++i) {
    for (let k = 0; k < n; ++k) {
      for (let c = 0; c < 3; ++c) {
        // a coefficient not finite (a file's NaN) left out: the scale every splat's coefficients are read by
        const v = Math.abs(getCoefficient(i, k, c));
        if (v > scale && Number.isFinite(v)) scale = v;
      }
    }
  }
  const data = new Uint8Array(splatCount * n * 4);
  const f = scale > 0 ? 127.5 / scale : 0;
  for (let i = 0; i < splatCount; ++i) {
    for (let k = 0; k < n; ++k) {
      const o = (i * n + k) * 4;
      for (let c = 0; c < 3; ++c) {
        const v = getCoefficient(i, k, c);
        data[o + c] = Number.isFinite(v) ? Math.round(v * f + 127.5) : 128;
      }
    }
  }
  return { degree, scale, data };
}

/** Coefficient `basis` (past the zeroth) of channel `channel` of splat `i` of `sh` */
export function getSplatShCoefficient(sh: SplatSh, i: number, basis: number, channel: number) {
  const n = ShBasisCounts[sh.degree];
  return (sh.data[(i * n + basis) * 4 + channel] / 255) * 2 * sh.scale - sh.scale;
}
