/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import type { TextureImage } from '@molstar/graphics/gl/renderable/util';
import { clamp, spline } from '@molstar/core/math/interpolate';
import { ValueCell } from '@molstar/core/util';
import type { Vec2 } from '@molstar/core/math/linear-algebra';

export interface ControlPoint {
  x: number;
  alpha: number;
}

export function getControlPointsFromString(s: string): ControlPoint[] {
  return s.split(/\s*,\s*/).map((p) => {
    const ps = p.split(/\s*:\s*/);
    return { x: parseFloat(ps[0]), alpha: parseFloat(ps[1]) };
  });
}

export function getControlPointsFromVec2Array(array: Vec2[]): ControlPoint[] {
  return array.map((v) => ({ x: v[0], alpha: v[1] }));
}

const TF_WIDTH = 256;

/** Build the 1D alpha curve (length 256) from sorted control points using a Catmull-Rom spline. */
function build1DAlphaArray(controlPoints: ControlPoint[], out?: Uint8Array): Uint8Array {
  const cp = [{ x: 0, alpha: 0 }, { x: 0, alpha: 0 }, ...controlPoints, { x: 1, alpha: 0 }, { x: 1, alpha: 0 }];

  const array = out ?? new Uint8Array(TF_WIDTH);
  if (out) array.fill(0);

  let k = 0;
  let x1: number, x2: number;
  let a0: number, a1: number, a2: number, a3: number;

  const il = controlPoints.length + 1;
  for (let i = 0; i < il; ++i) {
    x1 = cp[i + 1].x;
    x2 = cp[i + 2].x;

    a0 = cp[i].alpha;
    a1 = cp[i + 1].alpha;
    a2 = cp[i + 2].alpha;
    a3 = cp[i + 3].alpha;

    const jl = Math.round((x2 - x1) * TF_WIDTH);
    for (let j = 0; j < jl; ++j) {
      const t = j / jl;
      if (k < TF_WIDTH) {
        // the spline overshoots next to an opacity of 1: clamped, as the bytes would wrap around (256 to 0)
        array[k] = clamp(spline(a0, a1, a2, a3, t, 0.5) * 255, 0, 255);
        ++k;
      }
    }
  }
  return array;
}

export function createTransferFunctionTexture(
  controlPoints: ControlPoint[],
  texture?: ValueCell<TextureImage<Uint8Array>>,
): ValueCell<TextureImage<Uint8Array>> {
  // If updating an existing 2D texture, fall back to creating a fresh image so the size matches.
  const reuse = texture && texture.ref.value.height === 1 && texture.ref.value.array.length === TF_WIDTH;
  const array = reuse ? texture!.ref.value.array : new Uint8Array(TF_WIDTH);
  build1DAlphaArray(controlPoints, array);

  const textureImage = { array, width: TF_WIDTH, height: 1 };
  if (texture) {
    ValueCell.update(texture, textureImage);
    return texture;
  } else {
    return ValueCell.create(textureImage);
  }
}

/** The highest opacity per cell the pre-integration takes (its optical depth is finite) */
const MaxPreIntegratedAlpha = 0.999;

/**
 * The pre-integrated transfer function (`preIntegration`): for a ray step from value `i` to value `j` (rows `j`,
 * columns `i`, of 256 each, as the 1D table's bytes), the opacity of crossing one cell with the value varying
 * linearly between them, as an opacity per cell: the mean optical depth `-ln(1 - alpha)` over the values between,
 * `1 - exp(-mean)`. Independent of the step's length (the shader scales a cell's opacity to its step), so built once
 * per transfer function; on the diagonal the 1D table's own opacity. Narrow peaks crossed between two samples count,
 * which with the 1D table only show where a sample lands on them (slabs, and banding at few steps per cell)
 */
export function createPreIntegratedTransferTexture(
  alpha: Uint8Array,
  texture?: ValueCell<TextureImage<Uint8Array>>,
): ValueCell<TextureImage<Uint8Array>> {
  const n = TF_WIDTH;
  const depth = new Float64Array(n);
  for (let i = 0; i < n; ++i) depth[i] = -Math.log(1 - Math.min(alpha[i] / 255, MaxPreIntegratedAlpha));
  // the integral of the depth from value 0 to each value, linear between the table's values (trapezoids)
  const integral = new Float64Array(n);
  for (let i = 1; i < n; ++i) integral[i] = integral[i - 1] + (depth[i - 1] + depth[i]) / 2;
  const reuse = texture && texture.ref.value.array.length === n * n;
  const array = reuse ? texture!.ref.value.array : new Uint8Array(n * n);
  for (let j = 0; j < n; ++j) {
    for (let i = 0; i < n; ++i) {
      const mean = i === j ? depth[i] : (integral[j] - integral[i]) / (j - i);
      array[j * n + i] = Math.round((1 - Math.exp(-mean)) * 255);
    }
  }
  const image = { array, width: n, height: n };
  if (texture) {
    ValueCell.update(texture, image);
    return texture;
  }
  return ValueCell.create(image);
}
