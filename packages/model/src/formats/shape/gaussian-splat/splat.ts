/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type SplatFile, SplatRow, SplatRowBytes } from '@molstar/io/reader/splat/schema';
import {
  type GaussianSplat,
  type GaussianSplatData,
  type GaussianSplatReadOptions,
  setCenterRotation,
} from './data.js';

/**
 * The gaussian splats of antimatter15's `.splat`: linear scales (their logs here), the base color and opacity as bytes,
 * the rotation `w, x, y, z` as bytes `q * 128 + 128`; the frame of the 3DGS PLY it is made from (y-down). No
 * view-dependent color
 */
export function gaussianSplatsFromSplat(file: SplatFile, _options?: GaussianSplatReadOptions): GaussianSplatData {
  const { data, view } = file;
  return {
    count: file.count,
    upright: false,
    getSplat: (i: number, out: GaussianSplat) => {
      const o = i * SplatRowBytes;
      const p = o + SplatRow.position;
      const q = o + SplatRow.rotation;
      const rot = (k: number) => (data[q + k] - 128) / 128;
      setCenterRotation(
        out,
        view.getFloat32(p, true),
        view.getFloat32(p + 4, true),
        view.getFloat32(p + 8, true),
        rot(1),
        rot(2),
        rot(3),
        rot(0),
      );
      const s = o + SplatRow.scale;
      out.logScaleX = Math.log(view.getFloat32(s, true));
      out.logScaleY = Math.log(view.getFloat32(s + 4, true));
      out.logScaleZ = Math.log(view.getFloat32(s + 8, true));
      const c = o + SplatRow.color;
      out.r = data[c] / 255;
      out.g = data[c + 1] / 255;
      out.b = data[c + 2] / 255;
      out.opacity = data[c + 3] / 255;
    },
  };
}
