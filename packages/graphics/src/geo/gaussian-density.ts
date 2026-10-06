import type { Box3D, PositionData } from '@molstar/core/math/geometry';
import type { GaussianDensityProps } from '@molstar/core/math/geometry/gaussian-density';
import type { Mat4, Vec3, Vec2 } from '@molstar/core/math/linear-algebra';
import type { WebGLContext } from '@molstar/graphics/gl/webgl/context';
import type { Texture } from '@molstar/graphics/gl/webgl/texture';
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
  webgl: WebGLContext,
  texture?: Texture,
) {
  return _computeGaussianDensityTexture(webgl.isWebGL2 ? '3d' : '2d', position, box, radius, props, webgl, texture);
}

export function computeGaussianDensityTexture2d(
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
  webgl: WebGLContext,
  texture?: Texture,
) {
  return _computeGaussianDensityTexture('2d', position, box, radius, props, webgl, texture);
}

export function computeGaussianDensityTexture3d(
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
  webgl: WebGLContext,
  texture?: Texture,
) {
  return _computeGaussianDensityTexture('2d', position, box, radius, props, webgl, texture);
}

function _computeGaussianDensityTexture(
  type: '2d' | '3d',
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
  webgl: WebGLContext,
  texture?: Texture,
) {
  return Task.create('Gaussian Density', async (ctx) => {
    return type === '2d'
      ? GaussianDensityTexture2d(webgl, position, box, radius, false, props, texture)
      : GaussianDensityTexture3d(webgl, position, box, radius, props, texture);
  });
}
