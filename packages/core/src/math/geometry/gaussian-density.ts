/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Box3D, DensityData } from '../geometry.js';
import type { PositionData } from './common.js';
import { Task } from '@molstar/core/task/task';
import { GaussianDensityCPU } from './gaussian-density/cpu.js';

export const DefaultGaussianDensityProps = {
  resolution: 1,
  radiusOffset: 0,
  smoothness: 1.5,
};
export type GaussianDensityProps = typeof DefaultGaussianDensityProps;

export type GaussianDensityData = {
  radiusFactor: number;
} & DensityData;

export function computeGaussianDensity(
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
) {
  return Task.create('Gaussian Density', async (ctx) => {
    return await GaussianDensityCPU(ctx, position, box, radius, props);
  });
}
