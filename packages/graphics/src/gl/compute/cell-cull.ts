/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { KernelCode, StorageBufferSpec } from '../compute.js';
import type { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { TextureSpec, UniformSpec } from '../renderable/schema.js';
import type { Texture } from '../texture.js';

/**
 * GPU culling of merged render objects by the cells of their instance grids (`merged-render-objects.md` GPU culling):
 * `cell-cull` tests each cell (a thread per cell) and writes the levels it is drawn at; `instance-copy` copies each
 * instance of a drawn cell (a thread per instance) into the compacted instance buffers, a slot per level in its
 * segment's region (an atomic counter per level and segment); `cull-args` writes the instance counts into the
 * indirect records (their other words written once per layout, `merged-culler.ts`). Counters and records are
 * level-major: level `j` of segment `s` is `j * segmentCount + s`, so a level's records are one list; its region
 * starts at instance `j * instanceCount + instanceBase(s)`
 */
export const CellCullKernel = KernelCode('cell-cull', [64]);
export const InstanceCopyKernel = KernelCode('instance-copy', [64]);
export const CullArgsKernel = KernelCode('cull-args', [64]);

/** Levels the kernels handle, as `aSegmentLod` (4 components) */
export const MaxCullLevels = 4;

/**
 * A depth pyramid to test occlusion against on the GPU (`HiZPass`): its levels from `minLevel` side by side in one
 * `r32f` texture (level `i` of size `2^(levelCount - 1 - i)`, at the sum of the sizes before it), the camera it was built
 * with; depth reversed (the near plane at 1)
 */
export interface OcclusionPyramid {
  readonly texture: Texture;
  readonly levelCount: number;
  readonly minLevel: number;
  /** the viewport size the levels cover */
  readonly width: number;
  readonly height: number;
  readonly view: Mat4;
  readonly projection: Mat4;
  readonly near: number;
  readonly far: number;
}

export const CellCullSchema = {
  /** the frustum's planes (`[nx, ny, nz, constant]`) */
  uPlane0: UniformSpec('v4'),
  uPlane1: UniformSpec('v4'),
  uPlane2: UniformSpec('v4'),
  uPlane3: UniformSpec('v4'),
  uPlane4: UniformSpec('v4'),
  uPlane5: UniformSpec('v4'),
  uCameraPlane: UniformSpec('v4'),
  /** the object's distance range (`uLod` min and max) in x and y, z > 0 where it has one */
  uLodRange: UniformSpec('v4'),
  /** per level (column) its distance range in x and y */
  uLevels: UniformSpec('m4'),
  /** 0 without levels (all drawn as level 0) */
  uLevelCount: UniformSpec('i'),
  uCellCount: UniformSpec('i'),
  /** occlusion against the pyramid (`OcclusionPyramid`), its camera and layout */
  uHiZ: UniformSpec('b'),
  uHiZView: UniformSpec('m4'),
  uHiZProjection: UniformSpec('m4'),
  /** near, far */
  uHiZDepthRange: UniformSpec('v2'),
  uHiZViewport: UniformSpec('v2'),
  uHiZLevelCount: UniformSpec('i'),
  uHiZMinLevel: UniformSpec('i'),
  tHiZ: TextureSpec('texture', 'r', 'float', 'nearest'),
  /** per cell its sphere */
  cellSpheres: StorageBufferSpec('vec4<f32>', 'read'),
  /** per cell the levels it is drawn at, a bit each (0: culled) */
  cellLevels: StorageBufferSpec('u32', 'read-write'),
};

export const InstanceCopySchema = {
  uInstanceCount: UniformSpec('i'),
  uSegmentCount: UniformSpec('i'),
  /** the instances' data (`packInstanceData`) */
  instanceData: StorageBufferSpec('vec4<f32>', 'read'),
  /** per instance its cell and segment */
  instanceCells: StorageBufferSpec('vec2<u32>', 'read'),
  /** per segment its first instance */
  segmentBases: StorageBufferSpec('u32', 'read'),
  cellLevels: StorageBufferSpec('u32', 'read'),
  /** per (level, segment) the instances drawn, zeroed before */
  counters: StorageBufferSpec('atomic<u32>', 'read-write'),
  /** the compacted instanced attributes: `aTransform` (4 vec4), `aInstance`, `aSegment` (3), the sphere, level scales */
  transforms: StorageBufferSpec('vec4<f32>', 'read-write'),
  instances: StorageBufferSpec('f32', 'read-write'),
  segments: StorageBufferSpec('f32', 'read-write'),
  spheres: StorageBufferSpec('vec4<f32>', 'read-write'),
  lods: StorageBufferSpec('vec4<f32>', 'read-write'),
};

export const CullArgsSchema = {
  /** records: levels times segments */
  uRecordCount: UniformSpec('i'),
  /** words per record: 5 indexed (`drawIndexedIndirect`), 4 not (`drawIndirect`); the instance count is the second */
  uRecordWords: UniformSpec('i'),
  counters: StorageBufferSpec('u32', 'read'),
  args: StorageBufferSpec('u32', 'read-write'),
};
