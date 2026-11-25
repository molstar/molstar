/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SpzColorScale, type SpzFile, SpzShDims } from '@molstar/io/reader/spz/schema';
import {
  createShAccessor,
  type GaussianSplat,
  type GaussianSplatData,
  type GaussianSplatReadOptions,
  SH_C0,
  setCenterRotation,
} from './data.js';

function halfToFloat(h: number) {
  const s = h & 0x8000 ? -1 : 1;
  const e = (h >> 10) & 0x1f;
  const f = h & 0x3ff;
  if (e === 0) return s * 2 ** -14 * (f / 1024);
  if (e === 31) return f ? Number.NaN : s * Number.POSITIVE_INFINITY;
  return s * 2 ** (e - 15) * (1 + f / 1024);
}

/**
 * The gaussian splats of Niantic's `.spz` (versions 1 to 3, see `SpzFile`): stored right-up-back (upright), the PLY's
 * frame turned 180° about x. View-dependent color of degree 4 read to 3
 * (the coefficients ordered by degree)
 */
export function gaussianSplatsFromSpz(file: SpzFile, options: GaussianSplatReadOptions): GaussianSplatData {
  const { version, positions, alphas, colors, scales, rotations } = file;
  const fixed = 1 / (1 << file.fractionalBits);
  const position =
    version === 1
      ? (i: number, c: number) => {
          const o = (i * 3 + c) * 2;
          return halfToFloat(positions[o] | (positions[o + 1] << 8));
        }
      : (i: number, c: number) => {
          const o = (i * 3 + c) * 3;
          // 24-bit signed: sign-extended from bit 23
          return (((positions[o] | (positions[o + 1] << 8) | (positions[o + 2] << 16)) << 8) >> 8) * fixed;
        };
  const q = [0, 0, 0, 0];
  const rotation =
    version >= 3
      ? (i: number) => {
          const o = i * 4;
          let comp =
            (rotations[o] | (rotations[o + 1] << 8) | (rotations[o + 2] << 16) | (rotations[o + 3] << 24)) >>> 0;
          // the smallest three of 10 bits, the last index in the low bits
          const largest = comp >>> 30;
          let sum = 0;
          for (let k = 3; k >= 0; --k) {
            if (k === largest) continue;
            const mag = comp & 511;
            const neg = (comp >>> 9) & 1;
            comp >>>= 10;
            q[k] = ((Math.SQRT1_2 * mag) / 511) * (neg ? -1 : 1);
            sum += q[k] * q[k];
          }
          q[largest] = Math.sqrt(Math.max(0, 1 - sum));
        }
      : (i: number) => {
          const o = i * 3;
          q[0] = rotations[o] / 127.5 - 1;
          q[1] = rotations[o + 1] / 127.5 - 1;
          q[2] = rotations[o + 2] / 127.5 - 1;
          q[3] = Math.sqrt(Math.max(0, 1 - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]));
        };
  const color = (i: number, c: number) => 0.5 + SH_C0 * ((colors[i * 3 + c] / 255 - 0.5) / SpzColorScale);

  const dim = SpzShDims[file.shDegree];
  const { sh } = file;
  return {
    count: file.count,
    upright: true,
    getSplat: (i: number, out: GaussianSplat) => {
      rotation(i);
      setCenterRotation(out, position(i, 0), position(i, 1), position(i, 2), q[0], q[1], q[2], q[3]);
      out.logScaleX = scales[i * 3] / 16 - 10;
      out.logScaleY = scales[i * 3 + 1] / 16 - 10;
      out.logScaleZ = scales[i * 3 + 2] / 16 - 10;
      out.opacity = alphas[i] / 255;
      out.r = color(i, 0);
      out.g = color(i, 1);
      out.b = color(i, 2);
    },
    sh: createShAccessor(
      Math.min(options.shDegree, file.shDegree),
      (i, k, c) => (sh[(i * dim + k) * 3 + c] - 128) / 128,
    ),
  };
}
