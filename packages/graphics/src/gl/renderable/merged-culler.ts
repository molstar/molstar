/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import type { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import type { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import type { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import type { ComputeContext, ComputeKernel, KernelRequest, StorageBuffer } from '../compute.js';
import {
  CellCullKernel,
  CellCullSchema,
  CullArgsKernel,
  CullArgsSchema,
  InstanceCopyKernel,
  InstanceCopySchema,
  MaxCullLevels,
  type OcclusionPyramid,
} from '../compute/cell-cull.js';

/** A 1×1 pyramid texture, bound while the cells are not tested for occlusion */
function createEmptyHiZ(resources: GraphicsResources): Texture {
  const texture = resources.texture('image-float32', 'r', 'float', 'nearest');
  // loaded: storage of its own (a defined texture gets it lazily)
  texture.load({ array: new Float32Array(1), width: 1, height: 1 });
  return texture;
}
import { IndirectRecordBytes, type MultiDrawBaseData } from '../render-item.js';
import type { GraphicsResources } from '../context.js';
import type { Texture } from '../texture.js';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { type LodLevelsValue, LodView } from '../renderable.js';
import { InstanceDataVec4s, type Merged, type MergeableValues, packInstanceData } from './merged.js';

/**
 * GPU culling of a merged render object by the cells of its members' instance grids (`merged-render-objects.md` GPU
 * culling): per frame (`dispatch`, before the passes) the cells tested and the instances of the drawn ones copied into
 * compacted instanced attributes (`attributes`, bound in place of the merged object's), a region per level and
 * segment, drawn by an indirect record per level and segment (`getLists`); the programs stay the same. Members without
 * a grid are a cell of all their instances. The compacted attributes end with all instances in their own order
 * (`identityBase`): views not culled on the GPU (a scaled model, another view of the frame) draw the CPU culling's
 * lists shifted to it
 */
export interface MergedCuller {
  /** the compacted instanced attributes, by name */
  readonly attributes: {
    readonly aTransform: ValueCell<StorageBuffer>;
    readonly aInstance: ValueCell<StorageBuffer>;
    readonly aSegment: ValueCell<StorageBuffer>;
    readonly aSegmentSphere: ValueCell<StorageBuffer>;
    readonly aSegmentLod: ValueCell<StorageBuffer>;
  };
  /** the first instance of the identity region */
  readonly identityBase: number;
  /** after the merged values synced: the cells, records and instances written where they changed */
  sync: () => void;
  /**
   * culls with the planes, the levels' ranges for `lod`, and the occlusion of `pyramid` where given (the previous
   * frame's); false where the kernels
   * are not ready yet (compiling in the background)
   */
  dispatch: (cameraPlane: Plane3D, frustum: Frustum3D, lod: LodView, pyramid?: OcclusionPyramid) => boolean;
  /** per level its indirect list, with the uniforms of the level lists of the CPU culling (`levels`) */
  getLists: (levels: readonly MultiDrawBaseData[]) => MultiDrawBaseData[];
  destroy: () => void;
}

/** The culling's kernels, the same for every object: prepared once a culler is created, or by `Canvas3D.precompile` */
export const MergedCullKernels: readonly KernelRequest[] = [
  { code: CellCullKernel, schema: CellCullSchema },
  { code: InstanceCopyKernel, schema: InstanceCopySchema },
  { code: CullArgsKernel, schema: CullArgsSchema },
];

/** Floats per instance of the compacted attributes, by name */
const AttributeFloats = { aTransform: 16, aInstance: 1, aSegment: 3, aSegmentSphere: 4, aSegmentLod: 4 } as const;

const EmptyInt = new Int32Array(0);
const EmptyUint = new Uint32Array(0);

export function createMergedCuller(
  compute: ComputeContext,
  /** for a texture bound while there is no pyramid */
  resources: GraphicsResources,
  merged: Merged,
  members: readonly MergeableValues[],
  indexed: boolean,
): MergedCuller {
  const { values, segments } = merged;
  const stride = indexed ? IndirectRecordBytes.indexed : IndirectRecordBytes.vertex;
  const words = stride / 4;

  const createBuffer = (bytes: number, label: string) => compute.createBuffer(Math.max(16, bytes), label);
  const buffer = (label: string) => ValueCell.create(createBuffer(16, label));
  const instanceData = buffer('culled instance data');
  const instanceCells = buffer('culled instance cells');
  const segmentBases = buffer('culled segment bases');
  const cellSpheres = buffer('culled cell spheres');
  const cellLevels = buffer('culled cell levels');
  const counters = buffer('culled counters');
  const args = buffer('culled draw arguments');
  const attributes = {
    aTransform: buffer('culled transforms'),
    aInstance: buffer('culled instances'),
    aSegment: buffer('culled segments'),
    aSegmentSphere: buffer('culled segment spheres'),
    aSegmentLod: buffer('culled segment levels'),
  };

  const levelCountCell = ValueCell.create(0);
  const cellValues = {
    uPlane0: ValueCell.create(Vec4()),
    uPlane1: ValueCell.create(Vec4()),
    uPlane2: ValueCell.create(Vec4()),
    uPlane3: ValueCell.create(Vec4()),
    uPlane4: ValueCell.create(Vec4()),
    uPlane5: ValueCell.create(Vec4()),
    uCameraPlane: ValueCell.create(Vec4()),
    uLodRange: ValueCell.create(Vec4()),
    uLevels: ValueCell.create(Mat4()),
    uLevelCount: levelCountCell,
    uCellCount: ValueCell.create(0),
    uHiZ: ValueCell.create(false),
    uHiZView: ValueCell.create(Mat4()),
    uHiZProjection: ValueCell.create(Mat4()),
    uHiZDepthRange: ValueCell.create(Vec2()),
    uHiZViewport: ValueCell.create(Vec2()),
    uHiZLevelCount: ValueCell.create(0),
    uHiZMinLevel: ValueCell.create(0),
    tHiZ: ValueCell.create(createEmptyHiZ(resources)),
    cellSpheres,
    cellLevels,
  };
  const emptyHiZ = cellValues.tHiZ.ref.value;
  const copyValues = {
    uInstanceCount: ValueCell.create(0),
    uSegmentCount: ValueCell.create(0),
    instanceData,
    instanceCells,
    segmentBases,
    cellLevels,
    counters,
    transforms: attributes.aTransform,
    instances: attributes.aInstance,
    segments: attributes.aSegment,
    spheres: attributes.aSegmentSphere,
    lods: attributes.aSegmentLod,
  };
  const argsValues = {
    uRecordCount: ValueCell.create(0),
    uRecordWords: ValueCell.create(words),
    counters,
    args,
  };
  let kernels: ComputeKernel[] | undefined;
  let preparing = false;
  const getKernels = () => {
    if (kernels) return true;
    if (!compute.hasKernels(MergedCullKernels)) {
      if (!preparing) {
        preparing = true;
        compute.prepareKernels(MergedCullKernels).finally(() => (preparing = false));
      }
      return false;
    }
    kernels = [
      compute.createKernel(CellCullKernel, CellCullSchema, cellValues),
      compute.createKernel(InstanceCopyKernel, InstanceCopySchema, copyValues),
      compute.createKernel(CullArgsKernel, CullArgsSchema, argsValues),
    ];
    return true;
  };

  /** replaces the buffer of `cell` by one of at least `bytes` where it is smaller */
  const ensure = (cell: ValueCell<StorageBuffer>, bytes: number) => {
    if (cell.ref.value.byteLength >= bytes) return;
    cell.ref.value.destroy();
    ValueCell.update(cell, createBuffer(bytes, 'merged culling'));
  };

  let instanceCount = 0;
  let segmentCount = 0;
  let cellCount = 0;
  let levelCount = 0;
  let identityBase = 0;

  function syncLayout() {
    const lodLevels = members[0].lodLevels?.ref.value as LodLevelsValue | undefined;
    levelCount = lodLevels ? Math.min(lodLevels.length, MaxCullLevels) : 0;
    const levels = Math.max(1, levelCount);
    instanceCount = values.instanceCount.ref.value;
    segmentCount = segments.count;
    identityBase = levels * instanceCount;

    // the cells: per member those of its grid with instances, or one of all its instances without a grid
    const spheres: number[] = [];
    const cellOf = new Uint32Array(instanceCount * 2);
    for (let s = 0; s < segmentCount; ++s) {
      const member = members[s];
      const base = segments.instanceBases[s];
      const grid = member.instanceGrid.ref.value;
      if (grid.cellSize <= 1) {
        const b = member.boundingSphere.ref.value as Sphere3D;
        const c = spheres.length / 4;
        spheres.push(b.center[0], b.center[1], b.center[2], b.radius);
        for (let i = 0; i < segments.instanceCounts[s]; ++i) cellOf.set([c, s], (base + i) * 2);
        continue;
      }
      const { cellOffsets, cellSpheres: gridSpheres, cellCount: gridCells } = grid;
      for (let k = 0; k < gridCells; ++k) {
        if (cellOffsets[k + 1] === cellOffsets[k]) continue;
        const c = spheres.length / 4;
        for (let q = 0; q < 4; ++q) spheres.push(gridSpheres[k * 4 + q]);
        for (let i = cellOffsets[k]; i < cellOffsets[k + 1]; ++i) cellOf.set([c, s], (base + i) * 2);
      }
    }
    cellCount = spheres.length / 4;

    ensure(instanceData, instanceCount * InstanceDataVec4s * 16);
    ensure(instanceCells, instanceCount * 8);
    ensure(segmentBases, segmentCount * 4);
    ensure(cellSpheres, cellCount * 16);
    ensure(cellLevels, cellCount * 4);
    ensure(counters, levels * segmentCount * 4);
    ensure(args, levels * segmentCount * stride);
    for (const [name, floats] of Object.entries(AttributeFloats)) {
      ensure(attributes[name as keyof typeof AttributeFloats], (levels + 1) * instanceCount * floats * 4);
    }
    instanceCells.ref.value.write(cellOf);
    segmentBases.ref.value.write(Uint32Array.from(segments.instanceBases.slice(0, segmentCount)));
    cellSpheres.ref.value.write(new Float32Array(spheres));

    // the records' words but the instance count: the level's count, first vertex or index, the region's first instance
    const records = new Uint32Array(levels * segmentCount * words);
    for (let j = 0; j < levels; ++j) {
      for (let s = 0; s < segmentCount; ++s) {
        const member = members[s];
        const drawCount: number = member.drawCount.ref.value;
        // gridless members draw at full detail at each level they overlap, as the CPU culling (`cullSegment`)
        const gridless = member.instanceGrid.ref.value.cellSize <= 1;
        const memberLevels = member.lodLevels?.ref.value as LodLevelsValue | undefined;
        const o = (j * segmentCount + s) * words;
        records[o] = levelCount === 0 || gridless || !memberLevels ? drawCount : memberLevels[j][3];
        const first = j * instanceCount + segments.instanceBases[s];
        if (indexed) {
          records[o + 2] = segments.elementOffsets[s] / 4;
          records[o + 4] = first;
        } else {
          records[o + 2] = segments.drawBases[s];
          records[o + 3] = first;
        }
      }
    }
    args.ref.value.write(records);

    ValueCell.updateIfChanged(levelCountCell, levelCount);
    ValueCell.updateIfChanged(cellValues.uCellCount, cellCount);
    ValueCell.updateIfChanged(copyValues.uInstanceCount, instanceCount);
    ValueCell.updateIfChanged(copyValues.uSegmentCount, segmentCount);
    ValueCell.updateIfChanged(argsValues.uRecordCount, levels * segmentCount);
  }

  /** the levels' distance ranges for `lod` (`LodView.getRange`) into `levelRanges`, per level a column */
  const levelRanges = Mat4.zero();
  const range: number[] = [0, 0, 0];
  function setLevelRanges(lod: LodView) {
    Mat4.setZero(levelRanges);
    const lodLevels = members[0].lodLevels?.ref.value as LodLevelsValue | undefined;
    for (let j = 0; j < levelCount; ++j) {
      LodView.getRange(lodLevels!, j, lod, range);
      levelRanges[j * 4] = range[0];
      levelRanges[j * 4 + 1] = range[1];
    }
  }

  /** the instances' data for the kernel, and their own order after the regions (the identity region) */
  function syncInstances() {
    const data = new Float32Array(instanceCount * InstanceDataVec4s * 4);
    packInstanceData(values, data);
    instanceData.ref.value.write(data);
    for (const [name, floats] of Object.entries(AttributeFloats)) {
      const array = (values[name].ref.value as Float32Array).subarray(0, instanceCount * floats);
      attributes[name as keyof typeof AttributeFloats].ref.value.write(array, identityBase * floats * 4);
    }
  }

  /** what the cells and records are written from: the grids, the counts and the levels */
  const getLayoutKey = () => {
    let key = `${values.instanceCount.ref.value},${segments.count}`;
    for (const m of members) {
      key += `,${m.instanceGrid.ref.version},${m.boundingSphere.ref.version},${m.drawCount.ref.value},${m.lodLevels?.ref.version ?? -1}`;
    }
    return key;
  };
  const getInstanceVersions = () =>
    `${values.aTransform.ref.version},${values.aInstance.ref.version},${values.aSegment.ref.version},${values.aSegmentSphere.ref.version},${values.aSegmentLod.ref.version}`;
  let layoutKey = '';
  let instanceVersions = '';
  /** advanced where the layout or the instances changed: a dispatch culled with an older one is not current */
  let dataVersion = 0;
  function sync() {
    const key = getLayoutKey();
    const layoutChanged = key !== layoutKey;
    if (layoutChanged) {
      layoutKey = key;
      syncLayout();
    }
    const versions = getInstanceVersions();
    if (layoutChanged || versions !== instanceVersions) {
      instanceVersions = versions;
      syncInstances();
      dataVersion += 1;
    }
  }
  sync();
  // compiled from now on: a frame waiting for what is pending waits for them, culled on the CPU only meanwhile
  getKernels();

  const setPlane = (cell: ValueCell<Vec4>, plane: Plane3D) => {
    const v = cell.ref.value;
    Vec4.set(v, plane.normal[0], plane.normal[1], plane.normal[2], plane.constant);
    ValueCell.update(cell, v);
  };

  /**
   * What the last dispatch culled with: the planes, the distance range, the levels' ranges, `dataVersion` and the
   * pyramid's camera; while they are the same its compacted attributes and records are current (a still camera culls
   * nothing again, once tested against a pyramid of its own view)
   */
  const dispatched = new Float64Array(7 * 4 + 2).fill(NaN);
  const dispatchedLevels = Mat4.zero().fill(NaN);
  let dispatchedVersion = -1;
  const dispatchedPyramid = new Float64Array(32).fill(NaN);
  let dispatchedHiZ = false;
  function isPyramidDispatched(pyramid: OcclusionPyramid | undefined) {
    if (!pyramid) return !dispatchedHiZ;
    if (!dispatchedHiZ) return false;
    for (let i = 0; i < 16; ++i) {
      if (dispatchedPyramid[i] !== pyramid.view[i] || dispatchedPyramid[16 + i] !== pyramid.projection[i]) return false;
    }
    return true;
  }
  const getPlane = (i: number, cameraPlane: Plane3D, frustum: Frustum3D) => (i < 6 ? frustum[i as 0] : cameraPlane);
  function isDispatched(cameraPlane: Plane3D, frustum: Frustum3D, min: number, max: number) {
    if (dispatchedVersion !== dataVersion || dispatched[28] !== min || dispatched[29] !== max) return false;
    for (let i = 0; i < 16; ++i) if (dispatchedLevels[i] !== levelRanges[i]) return false;
    for (let i = 0; i < 7; ++i) {
      const { normal, constant } = getPlane(i, cameraPlane, frustum);
      const o = i * 4;
      if (
        dispatched[o] !== normal[0] ||
        dispatched[o + 1] !== normal[1] ||
        dispatched[o + 2] !== normal[2] ||
        dispatched[o + 3] !== constant
      ) {
        return false;
      }
    }
    return true;
  }
  function setDispatched(cameraPlane: Plane3D, frustum: Frustum3D, min: number, max: number) {
    for (let i = 0; i < 7; ++i) {
      const { normal, constant } = getPlane(i, cameraPlane, frustum);
      dispatched.set([normal[0], normal[1], normal[2], constant], i * 4);
    }
    dispatched[28] = min;
    dispatched[29] = max;
    Mat4.copy(dispatchedLevels, levelRanges);
    dispatchedVersion = dataVersion;
  }

  let lists: MultiDrawBaseData[] = [];
  let listsLevels: readonly MultiDrawBaseData[] | undefined;
  let listsKey = '';

  return {
    attributes,
    get identityBase() {
      return identityBase;
    },
    sync,
    dispatch: (cameraPlane, frustum, lod, pyramid) => {
      if (instanceCount === 0 || cellCount === 0) return false;
      const [min, max] = values.uLod.ref.value as Vec4;
      setLevelRanges(lod);
      if (isDispatched(cameraPlane, frustum, min, max) && isPyramidDispatched(pyramid)) return true;
      if (!getKernels()) return false;
      dispatchedHiZ = !!pyramid;
      ValueCell.updateIfChanged(cellValues.uHiZ, !!pyramid);
      if (pyramid) {
        dispatchedPyramid.set(pyramid.view, 0);
        dispatchedPyramid.set(pyramid.projection, 16);
        ValueCell.update(cellValues.uHiZView, Mat4.copy(cellValues.uHiZView.ref.value, pyramid.view));
        ValueCell.update(cellValues.uHiZProjection, Mat4.copy(cellValues.uHiZProjection.ref.value, pyramid.projection));
        ValueCell.update(
          cellValues.uHiZDepthRange,
          Vec2.set(cellValues.uHiZDepthRange.ref.value, pyramid.near, pyramid.far),
        );
        ValueCell.update(
          cellValues.uHiZViewport,
          Vec2.set(cellValues.uHiZViewport.ref.value, pyramid.width, pyramid.height),
        );
        ValueCell.updateIfChanged(cellValues.uHiZLevelCount, pyramid.levelCount);
        ValueCell.updateIfChanged(cellValues.uHiZMinLevel, pyramid.minLevel);
      }
      ValueCell.updateIfChanged(cellValues.tHiZ, pyramid?.texture ?? emptyHiZ);
      setPlane(cellValues.uPlane0, frustum[0]);
      setPlane(cellValues.uPlane1, frustum[1]);
      setPlane(cellValues.uPlane2, frustum[2]);
      setPlane(cellValues.uPlane3, frustum[3]);
      setPlane(cellValues.uPlane4, frustum[4]);
      setPlane(cellValues.uPlane5, frustum[5]);
      setPlane(cellValues.uCameraPlane, cameraPlane);
      ValueCell.update(cellValues.uLevels, Mat4.copy(cellValues.uLevels.ref.value, levelRanges));
      const hasLod = min !== 0 || max !== 0;
      ValueCell.update(cellValues.uLodRange, Vec4.set(cellValues.uLodRange.ref.value, min, max, hasLod ? 1 : 0, 0));

      const records = Math.max(1, levelCount) * segmentCount;
      const groups = (n: number) => Math.ceil(n / CellCullKernel.workgroupSize[0]);
      counters.ref.value.clear(0, records * 4);
      compute.beginComputePass({ label: 'merged cell culling' });
      kernels![0].dispatch(groups(cellCount));
      kernels![1].dispatch(groups(instanceCount));
      kernels![2].dispatch(groups(records));
      compute.endComputePass();
      setDispatched(cameraPlane, frustum, min, max);
      return true;
    },
    getLists: (levels) => {
      const key = `${segmentCount},${levelCount},${args.ref.value.id}`;
      if (levels === listsLevels && key === listsKey) return lists;
      listsLevels = levels;
      listsKey = key;
      lists = levels.slice(0, Math.max(1, levelCount)).map((level, j) => ({
        firsts: EmptyInt,
        counts: EmptyInt,
        offsets: EmptyInt,
        instanceCounts: EmptyInt,
        baseVertices: EmptyInt,
        baseInstances: EmptyUint,
        count: segmentCount,
        uniforms: level.uniforms,
        indirect: { buffer: args.ref.value, byteOffset: j * segmentCount * stride, indexed },
      }));
      return lists;
    },
    destroy: () => {
      if (kernels) for (const k of kernels) k.destroy();
      emptyHiZ.destroy();
      const cells = [instanceData, instanceCells, segmentBases, cellSpheres, cellLevels, counters, args];
      for (const cell of [...cells, ...Object.values(attributes)]) cell.ref.value.destroy();
    },
  };
}
