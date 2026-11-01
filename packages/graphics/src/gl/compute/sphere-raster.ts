/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createEmptyHiZ, createHiZValues, HiZSchema, type OcclusionPyramid, setHiZValues } from './cell-cull.js';
import type { GraphicsResources } from '../context.js';
import { Clip } from '@molstar/core/util/clip';
import {
  ClipStorageVec4s,
  getInstanceClipCount,
  getPixelClipCount,
  getStorageClipVersion,
  packStorageClipObjects,
} from '../clip-storage.js';
import { ValueCell } from '@molstar/core/util';
import type { Mutable } from '@molstar/core/util/type-helpers';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import {
  type ComputeContext,
  type ComputeKernel,
  KernelCode,
  type StorageBuffer,
  StorageBufferSpec,
} from '../compute.js';
import type { GraphicsRenderItem } from '../render-item.js';
import { type LodLevelsValue, LodView } from '../renderable.js';
import type { SpheresValues } from '../renderable/spheres.js';
import type { MergedRasterInput } from '../renderable/merged-culler.js';
import { UniformSpec } from '../renderable/schema.js';
import { packAffineRows } from '../affine-rows.js';
import { ReversedZeroToOne } from '../depth-range.js';
import { getScratchBuffer } from './scratch.js';

/**
 * Small opaque spheres drawn by kernels instead of quads (`sphere-compute-raster.md`): a pass appends the chunks of 64
 * spheres of each drawn instance at each level it is drawn at, `sphere-raster-cull` for unmerged objects (a thread
 * per instance: the frustum, the levels' ranges), `sphere-raster-regions` for merged ones culled on the GPU (a thread
 * per slot of the culler's compacted instances); `sphere-raster-args` writes their indirect dispatch, `sphere-raster`
 * (a workgroup per chunk, a thread per sphere) writes each pixel a sphere covers into the visibility buffer: per pixel
 * the nearest sphere's id and depth, one `atomicStoreMax` of `vec2u(id, depth)` (`atomicMinMax64`). A dispatch of
 * each per object: its data pages bound as they are (`RenderItem.getDataLayer`), its values packed into the kernels'
 * uniforms, its chunks in a region of the frame's (`SphereRasterTable`).
 */
export const SphereRasterCullKernel = KernelCode('sphere-raster-cull', [64]);
export const SphereRasterRegionsKernel = KernelCode('sphere-raster-regions', [64]);
export const SphereRasterArgsKernel = KernelCode('sphere-raster-args', [1]);
export const SphereRasterKernel = KernelCode('sphere-raster', [64]);

/** Levels of an object the kernels handle (as the merged culler, `MaxCullLevels`) */
export const MaxSphereRasterLevels = 4;

/** Bytes per pixel of the visibility buffer: the id and the depth (`vec2u`), 0 where empty */
export const SphereRasterVisibilityBytes = 8;

/** u32 words per chunk: (slot, first, end, level), (`aInstance`, group-instance offset, group offset, scale) */
export const SphereRasterChunkWords = 8;

/** What every pass reads */
const SphereRasterParams = {
  uModelView: UniformSpec('m4'),
  uProjection: UniformSpec('m4'),
  /** the camera plane in model space, as `uGlobal.uCameraPlane` */
  uCameraPlane: UniformSpec('v4'),
  uModelScale: UniformSpec('f'),
  /** the target's size in pixels */
  uViewport: UniformSpec('v2'),
  uNear: UniformSpec('f'),
  uInstanceCount: UniformSpec('i'),
  uVertexCount: UniformSpec('i'),
  /** per level (column) its range, fade and scale as `uLod`: min, max, overlap, scale (no fade or cull where 0) */
  uLevels: UniformSpec('m4'),
  /** per level its spheres (the first of `tPositionGroup`, in level order; unmerged objects) */
  uLevelCounts: UniformSpec('v4'),
  uLevelCount: UniformSpec('i'),
  /** the object's slot in the counters (`CounterWords` each) and the indirect arguments (8 words each) */
  uObjectIndex: UniformSpec('i'),
  /** the object's region of the frame's chunks */
  uChunkBase: UniformSpec('i'),
  uChunkCapacity: UniformSpec('i'),
};

export const SphereRasterCullSchema = {
  ...SphereRasterParams,
  /** the spheres' invariant bounding sphere, before the instance transforms */
  uBoundingSphere: UniformSpec('v4'),
  uGroupCount: UniformSpec('i'),
  /** occlusion of the instances against the view's pyramid (`SphereRasterView.pyramid`) */
  ...HiZSchema,
  /** the clip objects instances are clipped by (`getInstanceClipCount`), in `clipObjects` (`ClipStorageVec4s` each) */
  uClipObjectCount: UniformSpec('i'),
  clipObjects: StorageBufferSpec('vec4<f32>', 'read'),
  transforms: StorageBufferSpec('vec4<f32>', 'read'),
  instances: StorageBufferSpec('f32', 'read'),
  counter: StorageBufferSpec('atomic<u32>', 'read-write'),
  chunks: StorageBufferSpec('vec4<u32>', 'read-write'),
};

export const SphereRasterRegionsSchema = {
  ...SphereRasterParams,
  uSegmentCount: UniformSpec('i'),
  uRecordWords: UniformSpec('i'),
  /** the levels of the segments' scales (`tSegment`), 0 without levels */
  uSegmentLevels: UniformSpec('i'),
  /** `tSegment`'s first word in `segmentPage` */
  uSegmentOffset: UniformSpec('i'),
  /** the raster regions' first record and slot (`MergedRasterInput`) */
  uRecordBase: UniformSpec('i'),
  uSlotBase: UniformSpec('i'),
  records: StorageBufferSpec('u32', 'read'),
  instances: StorageBufferSpec('f32', 'read'),
  segments: StorageBufferSpec('vec4<f32>', 'read'),
  segmentPage: StorageBufferSpec('u32', 'read'),
  counter: StorageBufferSpec('atomic<u32>', 'read-write'),
  chunks: StorageBufferSpec('vec4<u32>', 'read-write'),
};

export const SphereRasterArgsSchema = {
  uObjectIndex: UniformSpec('i'),
  uChunkCapacity: UniformSpec('i'),
  /** 1: the chunks so far marked as drawn, before a phase appending more */
  uMark: UniformSpec('i'),
  counter: StorageBufferSpec('u32', 'read-write'),
  args: StorageBufferSpec('u32', 'read-write'),
};

export const SphereRasterSchema = {
  ...SphereRasterParams,
  /** the largest half size in pixels the raster draws, the rest left to the impostors (`rasterizesSphere`) */
  uMaxRadius: UniformSpec('f'),
  /** the rows of instances of one sphere (`appendChunks`), 64 a workgroup, not the chunks */
  uSingles: UniformSpec('i'),
  uPositionOffset: UniformSpec('i'),
  uSizeOffset: UniformSpec('i'),
  /** `uSizeIndexType`, -1 for a uniform size (`uSize`) */
  uSizeIndexType: UniformSpec('i'),
  uSize: UniformSpec('f'),
  uSizeFactor: UniformSpec('f'),
  /**
   * pixel clipping as the impostors' (`sphereImpostor`, `dClipSphere`): the clip objects (`getPixelClipCount`, in
   * `clipObjects`), whether they clip spheres whole by their centers instead (`dClipPrimitive`), whether a clipped
   * front shows the cap of a solid interior (`dSolidInterior`) or else the back (`uDoubleSided`); the model space of the
   * hits from the view's (`uModelView` inverted)
   */
  uClipObjectCount: UniformSpec('i'),
  uClipPrimitive: UniformSpec('i'),
  uSolidInterior: UniformSpec('i'),
  uDoubleSided: UniformSpec('i'),
  uInvModelView: UniformSpec('m4'),
  /**
   * a clipping layer (`dClipping`, `tClipping` from the data pages): its offset (-1 without), whether by instance
   * (`uInstanceGranularity`, else by group and instance); per sphere its flag exempting clip objects (`clipFlag`)
   */
  uClippingOffset: UniformSpec('i'),
  uClippingInstance: UniformSpec('i'),
  clipping: StorageBufferSpec('u32', 'read'),
  clipObjects: StorageBufferSpec('vec4<f32>', 'read'),
  transforms: StorageBufferSpec('vec4<f32>', 'read'),
  positions: StorageBufferSpec('u32', 'read'),
  sizes: StorageBufferSpec('u32', 'read'),
  counter: StorageBufferSpec('u32', 'read'),
  chunks: StorageBufferSpec('vec4<u32>', 'read'),
  vis: StorageBufferSpec('atomic<vec2<u32>>', 'read-write'),
};

/**
 * Pixels above which a sphere is drawn by its chunk's whole workgroup (its pixels strided over the threads), not by its
 * thread alone (`sphere-compute-raster.md` step 8: a thread per sphere loses to the impostors on large spheres)
 */
export const SphereRasterLargePixels = 64;

/** The kernels to prepare (`ComputeContext.prepareKernels`) */
export const SphereRasterKernels = [
  { code: SphereRasterCullKernel, schema: SphereRasterCullSchema },
  { code: SphereRasterRegionsKernel, schema: SphereRasterRegionsSchema },
  { code: SphereRasterArgsKernel, schema: SphereRasterArgsSchema },
  { code: SphereRasterKernel, schema: SphereRasterSchema },
];

/** The view the spheres are drawn in */
export interface SphereRasterView {
  /** the camera's view times the renderer's model matrix (`uModelView`) */
  readonly modelView: Mat4;
  /**
   * the camera's GL projection, perspective (orthographic views draw impostors); bound with reversed `[0, 1]` depth as
   * the WebGPU programs' (`ReversedZeroToOne`), so that a hit's depth is its `frag_depth`
   */
  readonly projection: Mat4;
  /** in model space (`uCameraPlane`) */
  readonly cameraPlane: Vec4;
  readonly modelScale: number;
  readonly near: number;
  /** the target's size in pixels */
  readonly width: number;
  readonly height: number;
  /** the largest half size on screen the raster draws, in pixels of the target */
  readonly maxRadius: number;
  /** the conversion of levels by size on screen (`Renderer`'s) */
  readonly lod: LodView;
  /**
   * the pyramid unmerged objects' instances are tested for occlusion against (that of GPU culling, the previous frame's,
   * `Renderer.setOcclusionPyramid`); merged objects' raster regions are of their culler's test
   */
  readonly pyramid?: OcclusionPyramid;
}

export interface SphereRasterObject {
  readonly item: GraphicsRenderItem;
  readonly values: SpheresValues;
  /** a merged object culled on the GPU this frame: its culler's compacted instances (`Renderable.getRasterInput`) */
  readonly culled?: MergedRasterInput;
}

/**
 * A frame's chunks: per object a region, the chunks two vec4u each (`SphereRasterChunkWords`); a sphere's id is its
 * chunk's index times 64 plus its lane, plus 1 (0 is no sphere), so ids are bounded by the frame's chunks
 */
export interface SphereRasterTable {
  readonly ranges: readonly SphereRasterRange[];
  readonly chunks: StorageBuffer;
}

/** An object's region of the frame's chunks, and what the resolve reads of it (`SphereResolveSchema`) */
export interface SphereRasterRange {
  readonly object: number;
  readonly chunkBase: number;
  readonly chunkCapacity: number;
  /** whether its spheres were drawn (not where its data does not fit or it has none) */
  readonly drawn: boolean;
  /** per level its range, fade and scale (`SphereRasterParams.uLevels`) */
  readonly levels: Mat4;
  /** the instances' transforms the chunks' slots index (affine rows) */
  readonly transforms: StorageBuffer | undefined;
}

/**
 * The object (its index in the drawn objects), instance slot, item, level and `aInstance` of `id`, by the chunks of
 * `table` as read back (`chunks`); undefined for none
 */
export function decodeSphereRasterId(table: SphereRasterTable, chunks: Uint32Array, id: number) {
  if (id === 0) return undefined;
  const chunk = Math.floor((id - 1) / 64);
  const range = table.ranges.find((r) => chunk >= r.chunkBase && chunk < r.chunkBase + r.chunkCapacity);
  if (!range) return undefined;
  const o = chunk * SphereRasterChunkWords;
  const item = chunks[o + 1] + ((id - 1) % 64);
  return { object: range.object, slot: chunks[o], item, level: chunks[o + 3], instance: chunks[o + 4] };
}

/** Whether the kernels can draw `values` at all: levels within `MaxSphereRasterLevels` */
export function canRasterSpheres(values: SpheresValues) {
  const levels = values.lodLevels.ref.value as LodLevelsValue | undefined;
  return (levels?.length ?? 0) <= MaxSphereRasterLevels;
}

type Cells<S> = { [k in keyof S]: ValueCell<any> };

type ObjectKernels = {
  /** the pass appending the chunks: `sphere-raster-cull` or `sphere-raster-regions` */
  chunk: ComputeKernel;
  args: ComputeKernel;
  raster: ComputeKernel;
  rasterSingles: ComputeKernel;
  chunkValues: Cells<typeof SphereRasterCullSchema> | Cells<typeof SphereRasterRegionsSchema>;
  rasterValues: Cells<typeof SphereRasterSchema>;
  argsValues: Cells<typeof SphereRasterArgsSchema>;
  /** unmerged objects: the transforms' rows and `aInstance`, packed where they changed */
  transforms?: StorageBuffer;
  instances?: StorageBuffer;
  instanceVersions: [number, number];
  /** the clip objects (`packStorageClipObjects`), written where they changed */
  clipObjects: StorageBuffer;
  clipVersion: string;
};

export interface SphereRaster {
  /**
   * Records the dispatches of `objects` into `vis` (`SphereRasterVisibilityBytes` per pixel of the view, cleared by the
   * caller) in a compute pass of its own; their render items updated (their data layers current), merged ones culled
   * for this frame
   */
  draw: (
    objects: readonly SphereRasterObject[],
    view: SphereRasterView,
    vis: StorageBuffer,
    /** a second phase (`second`) or `skipSecond` follows */
    twoPhase?: boolean,
  ) => SphereRasterTable;
  /**
   * The second phase of GPU culling (`MergedCuller.dispatchSecond`): the instances it found appended to the merged
   * objects' chunks and drawn into the same buffer, the chunks so far not again (the nearest stays, whichever phase)
   */
  second: () => void;
  /** No second phase after a `draw` announcing it */
  skipSecond: () => void;
  destroy: () => void;
}

/** An object's counters: chunk rows written, single rows written, those two drawn already, all rows asked for */
const CounterWords = 5;

/** The least chunks of a merged object to start with, 2 MB */
const InitialChunkCapacity = 65536;

export function createSphereRaster(compute: ComputeContext, resources: GraphicsResources): SphereRaster {
  const kernels = new Map<string, ObjectKernels>();
  /** bound while there is no pyramid */
  const emptyHiZ = createEmptyHiZ(resources);
  const range = [0, 0, 0];
  const levels = Mat4();
  const levelCounts = Vec4();
  let counter: StorageBuffer | undefined;
  let args: StorageBuffer | undefined;
  /** the objects of the last `draw` */
  const prepared: { k: ObjectKernels; threads: number; index: number; culled?: MergedRasterInput }[] = [];
  let readCounts: (() => void) | undefined;

  function getKernels(object: SphereRasterObject, vis: StorageBuffer, chunks: StorageBuffer) {
    const { item, culled } = object;
    const key = `${item.id}${culled ? ':merged' : ''}`;
    let k = kernels.get(key);
    if (k) return k;
    const common = {
      uModelView: ValueCell.create(Mat4()),
      uProjection: ValueCell.create(Mat4()),
      uCameraPlane: ValueCell.create(Vec4()),
      uModelScale: ValueCell.create(1),
      uViewport: ValueCell.create(Vec2()),
      uNear: ValueCell.create(0),
      uInstanceCount: ValueCell.create(0),
      uVertexCount: ValueCell.create(0),
      uLevels: ValueCell.create(Mat4()),
      uLevelCounts: ValueCell.create(Vec4()),
      uLevelCount: ValueCell.create(0),
      uObjectIndex: ValueCell.create(0),
      uChunkBase: ValueCell.create(0),
      uChunkCapacity: ValueCell.create(0),
    };
    const counterCell = ValueCell.create(counter!);
    const chunksCell = ValueCell.create(chunks);
    const transforms = culled ? undefined : compute.createBuffer(16 * 3, 'sphere raster transforms');
    const instances = culled ? undefined : compute.createBuffer(16, 'sphere raster instances');
    const clipObjects = compute.createBuffer(Clip.MaxObjectCount * ClipStorageVec4s * 16, 'sphere raster clip objects');
    const transformsCell = ValueCell.create(culled ? culled.transforms : transforms!);
    const chunkValues = culled
      ? {
          ...common,
          uSegmentCount: ValueCell.create(0),
          uRecordWords: ValueCell.create(0),
          uSegmentLevels: ValueCell.create(0),
          uSegmentOffset: ValueCell.create(0),
          uRecordBase: ValueCell.create(0),
          uSlotBase: ValueCell.create(0),
          records: ValueCell.create(culled.records),
          instances: ValueCell.create(culled.instances),
          segments: ValueCell.create(culled.segments),
          segmentPage: ValueCell.create(culled.records),
          counter: counterCell,
          chunks: chunksCell,
        }
      : {
          ...common,
          uBoundingSphere: ValueCell.create(Vec4()),
          uGroupCount: ValueCell.create(0),
          ...createHiZValues(emptyHiZ),
          uClipObjectCount: ValueCell.create(0),
          clipObjects: ValueCell.create(clipObjects),
          transforms: transformsCell,
          instances: ValueCell.create(instances!),
          counter: counterCell,
          chunks: chunksCell,
        };
    const rasterValues = {
      ...common,
      uMaxRadius: ValueCell.create(0),
      uSingles: ValueCell.create(0),
      uPositionOffset: ValueCell.create(0),
      uSizeOffset: ValueCell.create(0),
      uSizeIndexType: ValueCell.create(0),
      uSize: ValueCell.create(1),
      uSizeFactor: ValueCell.create(1),
      uClipObjectCount: ValueCell.create(0),
      uClipPrimitive: ValueCell.create(0),
      uSolidInterior: ValueCell.create(0),
      uDoubleSided: ValueCell.create(0),
      uInvModelView: ValueCell.create(Mat4()),
      uClippingOffset: ValueCell.create(-1),
      uClippingInstance: ValueCell.create(0),
      clipping: ValueCell.create(transforms ?? culled!.transforms),
      clipObjects: ValueCell.create(clipObjects),
      transforms: transformsCell,
      positions: ValueCell.create(transforms ?? culled!.transforms),
      sizes: ValueCell.create(transforms ?? culled!.transforms),
      counter: counterCell,
      chunks: chunksCell,
      vis: ValueCell.create(vis),
    };
    const argsValues = {
      uObjectIndex: common.uObjectIndex,
      uChunkCapacity: common.uChunkCapacity,
      uMark: ValueCell.create(0),
      counter: counterCell,
      args: ValueCell.create(args!),
    };
    k = {
      chunk: culled
        ? compute.createKernel(SphereRasterRegionsKernel, SphereRasterRegionsSchema, chunkValues as never)
        : compute.createKernel(SphereRasterCullKernel, SphereRasterCullSchema, chunkValues as never),
      args: compute.createKernel(SphereRasterArgsKernel, SphereRasterArgsSchema, argsValues),
      raster: compute.createKernel(SphereRasterKernel, SphereRasterSchema, rasterValues),
      rasterSingles: compute.createKernel(SphereRasterKernel, SphereRasterSchema, {
        ...rasterValues,
        uSingles: ValueCell.create(1),
      }),
      chunkValues,
      rasterValues,
      argsValues,
      transforms,
      instances,
      instanceVersions: [-1, -1],
      clipObjects,
      clipVersion: '',
    };
    kernels.set(key, k);
    return k;
  }

  /** unmerged objects: the instances' transforms (their rows) and `aInstance`, packed where they changed */
  function updateInstances(k: ObjectKernels, values: SpheresValues) {
    const count = values.instanceCount.ref.value;
    const { aTransform, aInstance } = values;
    if (k.instanceVersions[0] === aTransform.ref.version && k.instanceVersions[1] === aInstance.ref.version) return;
    const n = Math.max(1, count);
    if (k.transforms!.byteLength < n * 48) {
      k.transforms!.destroy();
      k.transforms = compute.createBuffer(n * 48, 'sphere raster transforms');
      ValueCell.update(k.rasterValues.transforms, k.transforms);
    }
    if (k.instances!.byteLength < Math.ceil(n / 4) * 16) {
      k.instances!.destroy();
      k.instances = compute.createBuffer(Math.ceil(n / 4) * 16, 'sphere raster instances');
      ValueCell.update(k.chunkValues.instances, k.instances);
    }
    if (count > 0) {
      const transform = aTransform.ref.value as Float32Array;
      k.transforms!.write(packAffineRows(transform, count, new Float32Array(count * 12)));
      k.instances!.write((aInstance.ref.value as Float32Array).subarray(0, count));
    }
    k.instanceVersions[0] = aTransform.ref.version;
    k.instanceVersions[1] = aInstance.ref.version;
  }

  /** the levels' uniforms into `levels` and `levelCounts` (the spheres of unmerged objects' levels), their count */
  function getLevels(lodLevels: LodLevelsValue | undefined, values: SpheresValues, view: SphereRasterView) {
    Mat4.setZero(levels);
    Vec4.set(levelCounts, 0, 0, 0, 0);
    if (!lodLevels || lodLevels.length === 0) {
      // the object's own range, as its draw block's `uLod`
      const lod = values.uLod.ref.value;
      levels[0] = lod[0];
      levels[1] = lod[1];
      levels[2] = lod[2];
      levels[3] = lod[3];
      levelCounts[0] = values.uVertexCount.ref.value;
      return 1;
    }
    const n = Math.min(lodLevels.length, MaxSphereRasterLevels);
    for (let j = 0; j < n; ++j) {
      LodView.getRange(lodLevels, j, view.lod, range);
      levels[j * 4] = range[0];
      levels[j * 4 + 1] = range[1];
      levels[j * 4 + 2] = range[2];
      levels[j * 4 + 3] = lodLevels[j][4];
      // vertices of the quads, 6 a sphere
      levelCounts[j] = lodLevels[j][3] / 6;
    }
    return n;
  }

  /**
   * The levels of an unmerged object: none without an instance grid, which is not culled and draws its spheres whole
   * (`createRenderable`'s `cull`)
   */
  function getUnmergedLevels(values: SpheresValues) {
    if (values.instanceGrid.ref.value.cellSize <= 1) return undefined;
    return values.lodLevels.ref.value as LodLevelsValue | undefined;
  }

  /**
   * Merged objects' chunks: every instance at its finest in view would be gigabytes, where views use a fraction. They
   * start small and follow the counts read back after a dispatch (`readChunkCounts`): a view needing more than half
   * doubles them for the next frame; one needing more than all loses spheres for that frame, the read being pending
   * work (the canvas draws again, an image pass renders again)
   */
  const chunkCapacities = new Map<string, number>();
  const chunkDispatches = new Map<string, number>();
  let chunkCountsRead = false;
  const getObjectKey = ({ item, culled }: SphereRasterObject) => `${item.id}${culled ? ':merged' : ''}`;
  function readChunkCounts(objects: readonly SphereRasterObject[], capacities: readonly number[]) {
    let dispatched = false;
    for (const o of objects) {
      if (!o.culled) continue;
      const key = getObjectKey(o);
      if (chunkDispatches.get(key) === o.culled.dispatchVersion) continue;
      chunkDispatches.set(key, o.culled.dispatchVersion);
      dispatched = true;
    }
    if (!dispatched || chunkCountsRead) return;
    chunkCountsRead = true;
    const read = counter!.readAsync().then((data) => {
      chunkCountsRead = false;
      const counts = new Uint32Array(data);
      for (let i = 0; i < objects.length; ++i) {
        const { culled } = objects[i];
        if (!culled) continue;
        const need = counts[i * CounterWords + 4];
        if (need * 2 <= capacities[i]) continue;
        chunkCapacities.set(getObjectKey(objects[i]), Math.min(need * 2, culled.chunkCapacity));
      }
    });
    resources.pending.track(read, 'sphere raster chunks');
  }

  /** the chunks of `object`: with each instance at each level it can be drawn at, merged objects' as they need */
  function getCapacity(object: SphereRasterObject) {
    const { values, culled } = object;
    if (culled) {
      const capacity =
        chunkCapacities.get(getObjectKey(object)) ?? Math.max(culled.instanceCount, InitialChunkCapacity);
      return Math.min(capacity, culled.chunkCapacity);
    }
    const lodLevels = getUnmergedLevels(values);
    let perInstance = 0;
    if (!lodLevels || lodLevels.length === 0) perInstance = Math.ceil(values.uVertexCount.ref.value / 64);
    else for (const l of lodLevels) perInstance += Math.ceil(l[3] / 6 / 64);
    return perInstance * values.instanceCount.ref.value;
  }

  return {
    draw: (objects, view, vis, twoPhase = false) => {
      const words = Math.max(1, objects.length) * CounterWords;
      if (!counter || counter.byteLength < words * 4) {
        counter?.destroy();
        args?.destroy();
        counter = compute.createBuffer(words * 4, 'sphere raster counters');
        args = compute.createBuffer(Math.max(1, objects.length) * 32, 'sphere raster args');
        for (const k of kernels.values()) {
          ValueCell.update(k.rasterValues.counter, counter);
          ValueCell.update(k.argsValues.args, args);
        }
      }
      counter.clear();

      // a region per object
      const ranges: Mutable<SphereRasterRange>[] = [];
      let chunkCount = 0;
      for (let i = 0; i < objects.length; ++i) {
        const chunkCapacity = getCapacity(objects[i]);
        ranges.push({
          object: i,
          chunkBase: chunkCount,
          chunkCapacity,
          drawn: false,
          levels: Mat4(),
          transforms: undefined,
        });
        chunkCount += chunkCapacity;
      }
      const chunks = getScratchBuffer(compute, 'sphere raster chunks', Math.max(1, chunkCount) * 32);

      prepared.length = 0;
      for (let i = 0; i < objects.length; ++i) {
        const object = objects[i];
        const { item, values, culled } = object;
        const instanceCount = culled ? culled.instanceCount : values.instanceCount.ref.value;
        if (instanceCount === 0 || values.uVertexCount.ref.value === 0 || ranges[i].chunkCapacity === 0) continue;
        const positions = item.getDataLayer?.('tPositionGroup');
        if (!positions) continue;
        const uniformSize = values.dSizeType.ref.value === 'uniform';
        const sizes = uniformSize ? positions : item.getDataLayer?.('tSize');
        if (!sizes) continue;
        const segmentLayer = culled ? item.getDataLayer?.('tSegment') : undefined;
        if (culled && !segmentLayer) continue;
        // a clipping layer's flags per sphere (pixel clipping): without its data the impostors draw the object
        const withClipping = getPixelClipCount(values) > 0 && values.dClipping.ref.value;
        const clipping = withClipping ? item.getDataLayer?.('tClipping') : undefined;
        if (withClipping && !clipping) continue;
        const k = getKernels(object, vis, chunks);
        if (!culled) updateInstances(k, values);

        const c = k.rasterValues;
        ValueCell.update(c.uModelView, Mat4.copy(c.uModelView.ref.value, view.modelView));
        ValueCell.update(c.uProjection, ReversedZeroToOne.projection(c.uProjection.ref.value, view.projection));
        ValueCell.update(c.uCameraPlane, Vec4.copy(c.uCameraPlane.ref.value, view.cameraPlane));
        ValueCell.updateIfChanged(c.uModelScale, view.modelScale);
        ValueCell.update(c.uViewport, Vec2.set(c.uViewport.ref.value, view.width, view.height));
        ValueCell.updateIfChanged(c.uNear, view.near);
        ValueCell.updateIfChanged(c.uInstanceCount, instanceCount);
        ValueCell.updateIfChanged(c.uVertexCount, values.uVertexCount.ref.value);
        const lodLevels = culled ? culled.lodLevels : getUnmergedLevels(values);
        const levelCount = getLevels(lodLevels, values, view);
        ValueCell.update(c.uLevels, Mat4.copy(c.uLevels.ref.value, levels));
        ValueCell.update(c.uLevelCounts, Vec4.copy(c.uLevelCounts.ref.value, levelCounts));
        ValueCell.updateIfChanged(c.uLevelCount, culled ? culled.levels : levelCount);
        ValueCell.updateIfChanged(c.uObjectIndex, i);
        ValueCell.updateIfChanged(c.uChunkBase, ranges[i].chunkBase);
        ValueCell.updateIfChanged(c.uChunkCapacity, ranges[i].chunkCapacity);
        // merged objects: their raster regions' instances whole (the culler's split per instance)
        ValueCell.updateIfChanged(c.uMaxRadius, culled ? 1e9 : view.maxRadius);
        ValueCell.updateIfChanged(c.uPositionOffset, positions.offset);
        ValueCell.updateIfChanged(c.uSizeOffset, uniformSize ? 0 : sizes.offset);
        ValueCell.updateIfChanged(c.uSizeIndexType, uniformSize ? -1 : values.uSizeIndexType.ref.value);
        ValueCell.updateIfChanged(c.uSize, values.uSize.ref.value);
        ValueCell.updateIfChanged(c.uSizeFactor, values.uSizeFactor.ref.value);
        ValueCell.updateIfChanged(c.positions, positions.buffer);
        ValueCell.updateIfChanged(c.sizes, sizes.buffer);
        ValueCell.updateIfChanged(c.chunks, chunks);
        ValueCell.updateIfChanged(c.vis, vis);
        if (culled) {
          const v = k.chunkValues as Cells<typeof SphereRasterRegionsSchema>;
          ValueCell.updateIfChanged(v.uSegmentCount, culled.segmentCount);
          ValueCell.updateIfChanged(v.uRecordWords, culled.recordWords);
          ValueCell.updateIfChanged(v.uSegmentLevels, culled.levelCount);
          ValueCell.updateIfChanged(v.uSegmentOffset, segmentLayer!.offset);
          ValueCell.updateIfChanged(v.uRecordBase, culled.recordBase);
          ValueCell.updateIfChanged(v.uSlotBase, culled.slotBase);
          ValueCell.updateIfChanged(v.records, culled.records);
          ValueCell.updateIfChanged(v.instances, culled.instances);
          ValueCell.updateIfChanged(v.segments, culled.segments);
          ValueCell.updateIfChanged(v.segmentPage, segmentLayer!.buffer);
          ValueCell.updateIfChanged(c.transforms, culled.transforms);
        } else {
          const v = k.chunkValues as Cells<typeof SphereRasterCullSchema>;
          ValueCell.update(
            v.uBoundingSphere,
            Vec4.copy(v.uBoundingSphere.ref.value, values.uInvariantBoundingSphere.ref.value),
          );
          ValueCell.updateIfChanged(v.uGroupCount, values.uGroupCount.ref.value);
          setHiZValues(v, view.pyramid, emptyHiZ);
          ValueCell.updateIfChanged(v.uClipObjectCount, getInstanceClipCount(values));
        }
        const clipVersion = getStorageClipVersion(values);
        if (k.clipVersion !== clipVersion) {
          k.clipVersion = clipVersion;
          const data = new Float32Array(Clip.MaxObjectCount * ClipStorageVec4s * 4);
          packStorageClipObjects(values, data);
          k.clipObjects.write(data);
        }
        ValueCell.updateIfChanged(c.uClipObjectCount, getPixelClipCount(values));
        ValueCell.updateIfChanged(c.uClipPrimitive, values.dClipPrimitive.ref.value ? 1 : 0);
        ValueCell.updateIfChanged(c.uSolidInterior, values.dSolidInterior.ref.value ? 1 : 0);
        ValueCell.updateIfChanged(c.uDoubleSided, values.uDoubleSided.ref.value ? 1 : 0);
        ValueCell.update(c.uInvModelView, Mat4.invert(c.uInvModelView.ref.value, view.modelView));
        ValueCell.updateIfChanged(c.uClippingOffset, clipping ? clipping.offset : -1);
        ValueCell.updateIfChanged(c.uClippingInstance, values.uInstanceGranularity.ref.value ? 1 : 0);
        ValueCell.updateIfChanged(c.clipping, clipping ? clipping.buffer : positions.buffer);
        const threads = culled ? culled.levels * instanceCount : instanceCount;
        prepared.push({ k, threads, index: i, culled });
        ranges[i].drawn = true;
        Mat4.copy(ranges[i].levels, levels);
        ranges[i].transforms = culled ? culled.transforms : k.transforms;
      }

      compute.beginComputePass({ label: 'sphere raster' });
      for (const { k, threads, index } of prepared) {
        const groups = Math.ceil(threads / 64);
        k.chunk.dispatch(Math.min(groups, 65535), Math.ceil(groups / 65535));
        k.args.dispatch(1);
        k.raster.dispatchIndirect(args!, index * 32);
        k.rasterSingles.dispatchIndirect(args!, index * 32 + 16);
      }
      compute.endComputePass();
      const drawnObjects = objects.slice();
      const capacities = ranges.map((r) => r.chunkCapacity);
      readCounts = () => readChunkCounts(drawnObjects, capacities);
      if (!twoPhase) readCounts();
      return { ranges, chunks };
    },
    second: () => {
      compute.beginComputePass({ label: 'sphere raster, second phase' });
      for (const { k, threads, index, culled } of prepared) {
        if (!culled?.secondDispatched) continue;
        ValueCell.updateIfChanged(k.argsValues.uMark, 1);
        k.args.dispatch(1);
        ValueCell.updateIfChanged(k.argsValues.uMark, 0);
        ValueCell.updateIfChanged(
          (k.chunkValues as Cells<typeof SphereRasterRegionsSchema>).records,
          culled.secondRecords,
        );
        const groups = Math.ceil(threads / 64);
        k.chunk.dispatch(Math.min(groups, 65535), Math.ceil(groups / 65535));
        k.args.dispatch(1);
        k.raster.dispatchIndirect(args!, index * 32);
        k.rasterSingles.dispatchIndirect(args!, index * 32 + 16);
      }
      compute.endComputePass();
      readCounts?.();
    },
    skipSecond: () => readCounts?.(),
    destroy: () => {
      for (const k of kernels.values()) {
        k.chunk.destroy();
        k.args.destroy();
        k.raster.destroy();
        k.rasterSingles.destroy();
        k.transforms?.destroy();
        k.instances?.destroy();
        k.clipObjects.destroy();
      }
      kernels.clear();
      emptyHiZ.destroy();
      counter?.destroy();
      args?.destroy();
    },
  };
}
