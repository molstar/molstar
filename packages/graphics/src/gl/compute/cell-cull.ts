/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { KernelCode, StorageBufferSpec } from '../compute.js';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { ValueCell } from '@molstar/core/util/value-cell';
import type { GraphicsResources } from '../context.js';
import { TextureSpec, UniformSpec } from '../renderable/schema.js';
import type { Texture } from '../texture.js';

/**
 * GPU culling of merged render objects by the cells of their instance grids (`merged-render-objects.md` GPU culling):
 * `cell-cull` tests each cell (a thread per cell) and writes the levels it is drawn at; `instance-copy` tests each
 * instance of a drawn cell (a thread per instance) the same way and copies it at the levels of both into the
 * compacted instance buffers, a slot per level in its
 * segment's region (an atomic counter per level and segment); `cull-args` writes the instance counts into the
 * indirect records (their other words written once per layout, `merged-culler.ts`). Counters and records are
 * level-major: level `j` of segment `s` is `j * segmentCount + s`, so a level's records are one list; its region
 * starts at instance `j * instanceCount + instanceBase(s)`
 */
export const CellCullKernel = KernelCode('cell-cull', [64]);
export const InstanceCopyKernel = KernelCode('instance-copy', [64]);
export const CullArgsKernel = KernelCode('cull-args', [64]);

/** Levels the kernels handle, as the segments' level scales (4 components, `tSegment`) */
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
  /**
   * advances with each build: a pyramid served again (not built in a frame, e.g. while work is pending) is the same
   * version, of the frame it was built in
   */
  readonly version: number;
}

/** Occlusion against a pyramid (`OcclusionPyramid`, `hiz_occlusion_wgsl`): its camera and layout */
export const HiZSchema = {
  uHiZ: UniformSpec('b'),
  uHiZView: UniformSpec('m4'),
  uHiZProjection: UniformSpec('m4'),
  /** near, far */
  uHiZDepthRange: UniformSpec('v2'),
  uHiZViewport: UniformSpec('v2'),
  uHiZLevelCount: UniformSpec('i'),
  uHiZMinLevel: UniformSpec('i'),
  tHiZ: TextureSpec('texture', 'r', 'float', 'nearest'),
};

type HiZValues = { [k in keyof typeof HiZSchema]: ValueCell<any> };

/** A 1×1 pyramid texture, bound while nothing is tested for occlusion */
export function createEmptyHiZ(resources: GraphicsResources): Texture {
  const texture = resources.texture('image-float32', 'r', 'float', 'nearest');
  // loaded: storage of its own (a defined texture gets it lazily)
  texture.load({ array: new Float32Array(1), width: 1, height: 1 });
  return texture;
}

/** The values of `HiZSchema` for `pyramid`; none: the test off, `empty` bound */
export function setHiZValues(values: HiZValues, pyramid: OcclusionPyramid | undefined, empty: Texture) {
  ValueCell.updateIfChanged(values.uHiZ, !!pyramid);
  if (pyramid) {
    ValueCell.update(values.uHiZView, Mat4.copy(values.uHiZView.ref.value, pyramid.view));
    ValueCell.update(values.uHiZProjection, Mat4.copy(values.uHiZProjection.ref.value, pyramid.projection));
    ValueCell.update(values.uHiZDepthRange, Vec2.set(values.uHiZDepthRange.ref.value, pyramid.near, pyramid.far));
    ValueCell.update(values.uHiZViewport, Vec2.set(values.uHiZViewport.ref.value, pyramid.width, pyramid.height));
    ValueCell.updateIfChanged(values.uHiZLevelCount, pyramid.levelCount);
    ValueCell.updateIfChanged(values.uHiZMinLevel, pyramid.minLevel);
  }
  ValueCell.updateIfChanged(values.tHiZ, pyramid?.texture ?? empty);
}

/** The cells of `HiZSchema`, the test off */
export function createHiZValues(empty: Texture): HiZValues {
  return {
    uHiZ: ValueCell.create(false),
    uHiZView: ValueCell.create(Mat4()),
    uHiZProjection: ValueCell.create(Mat4()),
    uHiZDepthRange: ValueCell.create(Vec2()),
    uHiZViewport: ValueCell.create(Vec2()),
    uHiZLevelCount: ValueCell.create(0),
    uHiZMinLevel: ValueCell.create(0),
    tHiZ: ValueCell.create(empty),
  };
}

/** What both kernels test a sphere against: the frustum, the object's and the levels' distance ranges, the pyramid */
const SphereTestSchema = {
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
  ...HiZSchema,
  /**
   * the clip objects of an object without a clipping layer (`getCullClipCount`), in `segmentData` from vec4
   * `uClipBase` (`clip_storage_wgsl`): spheres they clip whole culled (`clipSphereClass`, `clip-culling.md`)
   */
  uClipObjectCount: UniformSpec('i'),
  uClipBase: UniformSpec('i'),
  /**
   * the object's transient transform (`uObjectTransform`) and its largest scale: every sphere tested (a cell's, an
   * instance's) moved by it first, as the draws move the object
   */
  uObjectTransform: UniformSpec('m4'),
  uObjectScale: UniformSpec('f'),
  /**
   * the region of the ghost of what the clip objects remove (`Renderer.renderClipGhost`), 0 for none: the cells and
   * instances in view the clip objects reach copied into it too, at one level (the coarsest), in the first phase
   */
  uGhostRegion: UniformSpec('i'),
  /**
   * per segment three vec4s: its invariant sphere (an instance's sphere by its transform: tested as its cell's, drawn
   * at the levels of both); its first instance (the bits of a u32 in x) and its largest sphere's radius (y); its level
   * scales (`tSegment`'s); after the segments the clip objects (`uClipBase`)
   */
  segmentData: StorageBufferSpec('vec4<f32>', 'read'),
};

export const CellCullSchema = {
  ...SphereTestSchema,
  uCellCount: UniformSpec('i'),
  /** per cell its sphere */
  cellSpheres: StorageBufferSpec('vec4<f32>', 'read'),
  /** per cell the levels it is drawn at, a bit each (0: culled) */
  cellLevels: StorageBufferSpec('u32', 'read-write'),
};

export const InstanceCopySchema = {
  ...SphereTestSchema,
  uInstanceCount: UniformSpec('i'),
  uSegmentCount: UniformSpec('i'),
  /**
   * the sphere raster (`sphere-raster.ts`): its largest radius in pixels (0: no raster regions), the pixels per Å of
   * radius over the distance (the projection's focal scale), the near plane; its regions after the levels' (`uRasterLevels`
   * of them, the levels' count)
   */
  uRasterRadius: UniformSpec('f'),
  uRasterScale: UniformSpec('f'),
  uRasterNear: UniformSpec('f'),
  uRasterLevels: UniformSpec('i'),
  /** instance clipping (`getInstanceClipCount`): the instances' centers tested against the clip objects too */
  uClipInstance: UniformSpec('i'),
  /**
   * 1 in the second phase: the instances the first left out (the top bit of their cell, written by the first), appended
   * after its (the counters go on)
   */
  uPhase: UniformSpec('i'),
  /** per instance its cell, the top bit whether the first phase drew it, and segment */
  instanceCells: StorageBufferSpec('vec2<u32>', 'read-write'),
  cellLevels: StorageBufferSpec('u32', 'read'),
  /** per (level, segment) the instances drawn, zeroed before */
  counters: StorageBufferSpec('atomic<u32>', 'read-write'),
  /** the regions' instance indices, written (`sCulledIndices`; its identity region written once per layout) */
  indices: StorageBufferSpec('u32', 'read-write'),
  /** the instances' data, read by index (`CulledData`: the transforms' rows first) */
  data: StorageBufferSpec('vec4<f32>', 'read'),
};

export const CullArgsSchema = {
  /** records: levels times segments */
  uRecordCount: UniformSpec('i'),
  /** words per record: 5 indexed (`drawIndexedIndirect`), 4 not (`drawIndirect`); the instance count is the second */
  uRecordWords: UniformSpec('i'),
  /** 1 in the second phase: `args` holds the first's counts */
  uPhase: UniformSpec('i'),
  counters: StorageBufferSpec('u32', 'read'),
  args: StorageBufferSpec('u32', 'read-write'),
  /** the second phase's records, as `args` (its other words written once per layout) */
  second: StorageBufferSpec('u32', 'read-write'),
};
