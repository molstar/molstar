/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Mat4, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import type { GraphicsContext } from '../../context.js';
import type { Texture } from '../../texture.js';

/**
 * Marching cubes on a volume in a 2d atlas texture with fragment shaders (active voxels, histogram pyramid): an
 * implementation a backend registers (WebGL, `gl/webgl/compute/marching-cubes`). Backends with compute kernels use
 * `computeChunkedIsosurface`; neutral code chooses by capability, see `getAtlasIsosurfaceExtraction`.
 */

export interface AtlasIsosurfaceInput {
  /** the volume in a 2d atlas of slices, the value in `channel` */
  readonly texture: Texture;
  readonly channel: 'red' | 'alpha';
  readonly gridDim: Vec3;
  readonly gridTexDim: Vec3;
  readonly gridDataDim: Vec3;
  readonly gridTexScale: Vec2;
  readonly transform: Mat4;
  readonly isoValue: number;
  /** flips the winding and normals (negative iso values) */
  readonly invert: boolean;
  /** the group packed in the rgb of `texture`, else the cell index */
  readonly packedGroup: boolean;
  readonly axisOrder: Vec3;
  /** one group per triangle */
  readonly constantGroup: boolean;
}

export interface AtlasIsosurfaceOutput {
  readonly vertexTexture: Texture;
  readonly groupTexture: Texture;
  readonly normalTexture: Texture;
  readonly vertexCount: number;
}

export interface AtlasIsosurfaceExtraction {
  /** whether it runs on `graphics` (its backend) */
  isSupported(graphics: GraphicsContext): boolean;
  /** the isosurface, into `textures` when given */
  extract(
    graphics: GraphicsContext,
    input: AtlasIsosurfaceInput,
    textures?: { vertex: Texture; group: Texture; normal: Texture },
  ): AtlasIsosurfaceOutput;
}

const extractions: AtlasIsosurfaceExtraction[] = [];

/** Called by backends that implement it, when they create a context */
export function registerAtlasIsosurfaceExtraction(extraction: AtlasIsosurfaceExtraction) {
  if (!extractions.includes(extraction)) extractions.push(extraction);
}

/** The atlas isosurface extraction of the backend of `graphics`, if it has one */
export function getAtlasIsosurfaceExtraction(graphics: GraphicsContext): AtlasIsosurfaceExtraction | undefined {
  return extractions.find((e) => e.isSupported(graphics));
}
