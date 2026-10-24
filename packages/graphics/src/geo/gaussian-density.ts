/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Box3D, PositionData } from '@molstar/core/math/geometry';
import type { GaussianDensityProps } from '@molstar/core/math/geometry/gaussian-density';
import type { Mat4, Vec3, Vec2 } from '@molstar/core/math/linear-algebra';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Texture } from '@molstar/graphics/gl/texture';
import { Task } from '@molstar/core/task/task';
import { GaussianDensityTexture2d, GaussianDensityTexture3d } from './gaussian-density/gpu.js';

export type GaussianDensityTextureData = {
  radiusFactor: number;
  resolution: number;
  maxRadius: number;
  transform: Mat4;
  texture: Texture;
  bbox: Box3D;
  gridDim: Vec3;
  gridTexDim: Vec3;
  gridDataDim: Vec3;
  gridTexScale: Vec2;
};

export function computeGaussianDensityTexture(
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
  graphics: GraphicsContext,
  texture?: Texture,
) {
  return _computeGaussianDensityTexture(
    graphics.capabilities.texture3d ? '3d' : '2d',
    position,
    box,
    radius,
    props,
    graphics,
    texture,
  );
}

export function computeGaussianDensityTexture2d(
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
  graphics: GraphicsContext,
  texture?: Texture,
) {
  return _computeGaussianDensityTexture('2d', position, box, radius, props, graphics, texture);
}

export function computeGaussianDensityTexture3d(
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
  graphics: GraphicsContext,
  texture?: Texture,
) {
  return _computeGaussianDensityTexture('3d', position, box, radius, props, graphics, texture);
}

function _computeGaussianDensityTexture(
  type: '2d' | '3d',
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
  graphics: GraphicsContext,
  texture?: Texture,
) {
  return Task.create('Gaussian Density', async (ctx) => {
    return type === '2d'
      ? GaussianDensityTexture2d(graphics, position, box, radius, false, props, texture)
      : GaussianDensityTexture3d(graphics, position, box, radius, props, texture);
  });
}
