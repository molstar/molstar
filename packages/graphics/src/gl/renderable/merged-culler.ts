/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Clip } from '@molstar/core/util/clip';
import { ClipStorageVec4s, getCullClipCount, getInstanceClipCount, packStorageClipObjects } from '../clip-storage.js';
import { ValueCell } from '@molstar/core/util';
import { isDebugMode } from '@molstar/core/util/debug';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import type { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import type { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
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
  createEmptyHiZ,
  createHiZValues,
  setHiZValues,
} from '../compute/cell-cull.js';
import { getDispatchSize } from '../compute/prefix-sum.js';
import {
  createSceneCellBuilder,
  getSceneGrid,
  type SceneCellBuilder,
  SceneCellKernels,
  type SceneGrid,
} from '../compute/scene-cells.js';
import { type CulledPoseExpander, CulledPosesKernels, createCulledPoseExpander } from '../compute/culled-poses.js';
import { createU32RunWriter, type U32RunWriter, U32RunsKernels } from '../compute/u32-runs.js';
import { IndirectRecordBytes, type MultiDrawBaseData } from '../render-item.js';
import type { ItemIndexPattern } from '../quad-indices.js';
import type { GraphicsResources } from '../context.js';
import { type CullBatch, type LodLevelsValue, LodView } from '../renderable.js';
import { type CulledValues, type Merged, type MergeableValues, type MergedSegments } from './merged.js';
import { getObjectTransforms, getTransformParts, hasInstanceCells, type TransformValues } from '../transform-parts.js';
import { CulledData, getCulledSegmentData } from '../culled-data.js';
import { type DrawnBoundsValues, getDrawnBoundingSphere } from './util.js';
import type { InstanceGrid } from '@molstar/core/math/geometry/instance-grid';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';

/**
 * How often cullers split their cells (`layouts`), moved them as their grids moved (`moves`), built them again for
 * changed poses alone (`poses`) and dispatched a kernel in more than one row of workgroups (`rows`): for tests
 */
export const MergedCullerCounts = { layouts: 0, moves: 0, rows: 0, poses: 0 };
/**
 * The workgroups per dimension the kernels' dispatches go up to before a further row (`getDispatchSize`): WebGPU's
 * guaranteed limit; the instances `writeCulledData` writes at a time; the instances a
 * scene cell holds on average (`getSceneGrid`: the cells cost a thread each, an instance in a culled cell skips its
 * own test; every instance of a drawn cell is tested again, its sphere and the pyramid, so the cells' size decides the
 * work skipped, not what is drawn). Tests lower them to split small dispatches and writes, and to make many cells
 */
export const MergedCullerLimits = {
  maxWorkgroupsPerDimension: 65535,
  dataBlockInstances: 1 << 18,
  sceneCellInstances: 32,
};

/** `kernel` over `n` items, in rows of workgroups past the per-dimension limit (the kernels linearize the id) */
function dispatchItems(kernel: ComputeKernel, n: number) {
  const { x, y } = getDispatchSize(
    Math.ceil(n / CellCullKernel.workgroupSize[0]),
    MergedCullerLimits.maxWorkgroupsPerDimension,
  );
  if (y > 1) MergedCullerCounts.rows += 1;
  kernel.dispatch(x, y);
}
const IdentityPlacement = Mat4.identity();
const moveScratch = Mat4();

/**
 * GPU culling of a merged render object by the cells of its members' instance grids: per
 * frame (`dispatch`, before the passes) the cells tested and the indices of the drawn ones' instances
 * written into regions (`sCulledIndices`), a region per level and segment, drawn by an indirect record per level and
 * segment (`getLists`); the culled programs read each instance's data by its index (`CulledSchema`, `CulledData`).
 * Members without a grid are a cell of all their instances. The indices end with all instances in their own order
 * (`identityBase`): views not culled on the GPU (a scaled model, another view of the frame) draw the CPU culling's
 * lists shifted to it
 */
export interface MergedCuller {
  /**
   * a page store's culler (`merged-paged-culler.ts`): its instances composed on the GPU alone, views not culled there
   * drawing none
   */
  readonly paged?: boolean;
  /** what the culled programs read (`CulledSchema`), bound in place of the instanced attributes */
  readonly values: CulledValues;
  /** the first instance of the identity region */
  readonly identityBase: number;
  /**
   * The identity region's instances in `order` (slot `k` instance `order[k]`, all instances once), in place of their
   * own: the CPU culling's lists then ranges of that order (members' grids built after the data, in its order); none:
   * their own order again. Kept through layouts of the same instance count, written where the layout is current
   */
  setIdentityOrder: (order: Uint32Array | undefined) => void;
  /** after the merged values synced: the cells, records and instances written where they changed */
  sync: () => void;
  /** advanced by `sync` where the layout or the instances changed, and by a dispatch where the clip objects did */
  readonly dataVersion: number;
  /** whether the kernels are ready (compiled), their compile started where not (tracked as pending work) */
  areKernelsReady: () => boolean;
  /**
   * culls with the planes, the levels' ranges for `lod`, and the occlusion of `pyramid` where given (the previous
   * frame's); false where the kernels
   * are not ready yet (compiling in the background); with `batch` the clear and the dispatches recorded there
   */
  dispatch: (
    cameraPlane: Plane3D,
    frustum: Frustum3D,
    lod: LodView,
    pyramid?: OcclusionPyramid,
    raster?: MergedRasterParams,
    batch?: CullBatch,
    /** the ghost of what the clip objects remove too, into a region of its own (`getGhostLists`) */
    ghost?: boolean,
  ) => boolean;
  /** whether the last `dispatch` culled the ghost */
  readonly ghostDispatched: boolean;
  /**
   * the ghost's list (`Renderer.renderClipGhost`): its region at the coarsest level, at every distance (the uniforms of
   * the CPU culling's coarsest level list, `levels`, but its range)
   */
  getGhostLists: (levels: readonly MultiDrawBaseData[]) => MultiDrawBaseData[];
  /** whether `dispatchSecond` is to run for the last `dispatch`: it culled anew, against a pyramid */
  readonly secondPending: boolean;
  /**
   * The second phase of the last `dispatch` (no pass open): the instances it left out tested against `pyramid`, of the
   * frame's depth so far with the frame's camera, those visible for `getSecondLists`; what the first drew stays hidden
   * or not as it is, its depth part of the frame's, so what this leaves out is hidden in the finished frame. Afterwards
   * `getLists` draw both. No pyramid: nothing, the frame stays as the first phase left it. `batch` as `dispatch`'s
   */
  dispatchSecond: (pyramid: OcclusionPyramid | undefined, batch?: CullBatch) => boolean;
  /** whether `dispatchSecond` ran for the last `dispatch` */
  readonly secondDispatched: boolean;
  /** `getLists` for the instances of the second phase alone */
  getSecondLists: (levels: readonly MultiDrawBaseData[], withRaster?: boolean) => MultiDrawBaseData[];
  /** per level its indirect list, with the uniforms of the level lists of the CPU culling (`levels`) */
  getLists: (levels: readonly MultiDrawBaseData[], withRaster?: boolean) => MultiDrawBaseData[];
  /**
   * what the sphere raster draws from (`sphere-raster.ts`): the instances of the raster regions of the last `dispatch`;
   * undefined where it had no raster parameters
   */
  getRasterInput: () => MergedRasterInput | undefined;
  /** whether the lists draw only the records with instances (`compactRecords`) */
  readonly recordsCompact: boolean;
  /** advances where the records the lists draw changed: recordings of them no longer hold */
  readonly recordsVersion: number;
  destroy: () => void;
}

/**
 * The compacted instances of a dispatch for kernels drawing them (the sphere raster): a region per level and segment,
 * its count, vertex count and first vertex in `records` (the indirect records, level-major: `j * segmentCount + s`)
 */
export interface MergedRasterInput {
  /** the regions' instance indices (`sCulledIndices`) */
  readonly indices: StorageBuffer;
  /** the instances' data (`CulledData`): the transforms' affine rows first, three vec4s an instance */
  readonly data: StorageBuffer;
  readonly records: StorageBuffer;
  /** the second phase's records (`MergedCuller.dispatchSecond`), as `records` */
  readonly secondRecords: StorageBuffer;
  /** whether the second phase ran for this dispatch */
  readonly secondDispatched: boolean;
  /** u32 words a record (the first instance its last) */
  readonly recordWords: number;
  /** the levels' lists (at least 1), the levels of the first member (0 without) */
  readonly levels: number;
  readonly levelCount: number;
  readonly lodLevels: LodLevelsValue | undefined;
  readonly instanceCount: number;
  readonly segmentCount: number;
  /** chunks of 64 spheres of all regions full (each segment's instances at each level) */
  readonly chunkCapacity: number;
  /** advances with each dispatch: what the regions hold changed */
  readonly dispatchVersion: number;
  /** the raster regions' first record and first slot (after the levels' regions) */
  readonly recordBase: number;
  readonly slotBase: number;
}

/**
 * The sphere raster's split of a dispatch (the split per instance): instances whose largest
 * sphere is at most `radius` pixels at their nearest point, in front of the `near` plane, go into raster regions after
 * the levels' (`MergedRasterInput`); the impostors draw the rest (`getLists` without them)
 */
export interface MergedRasterParams {
  readonly radius: number;
  /** pixels per Å of radius over the distance: the projection's focal scale times half the viewport */
  readonly scale: number;
  readonly near: number;
}

/** The culling's kernels, the same for every object: prepared once a culler is created, or by `Canvas3D.precompile` */
export const MergedCullKernels: readonly KernelRequest[] = [
  { code: CellCullKernel, schema: CellCullSchema },
  { code: InstanceCopyKernel, schema: InstanceCopySchema },
  { code: CullArgsKernel, schema: CullArgsSchema },
  ...SceneCellKernels,
  ...CulledPosesKernels,
  ...U32RunsKernels,
];

/**
 * The transforms of `values` (`n` instances) by their encoding (`getTransformParts`: rows of matrices, or poses read
 * by index), `aInstance`, each instance's segment and each segment's groups (`segments`, as `aSegment` holds them per
 * instance) as `CulledData`, into `data` (grown where it is too small): written a block at a time, not composed into
 * one array first (millions of instances); the poses' slots expanded on the GPU (`expander`), and the words that are
 * runs (`CulledDataLayout.runs`: slot words, `aInstance`, segment indices) written there (`runWriter`), submitted before
 */
function writeCulledData(
  data: ValueCell<StorageBuffer>,
  ensure: (cell: ValueCell<StorageBuffer>, bytes: number) => void,
  values: TransformValues & { readonly aInstance: ValueCell<unknown> },
  segments: MergedSegments,
  n: number,
  expander: CulledPoseExpander,
  runWriter: U32RunWriter,
  /** the poses changed alone (`isPosesChange`): their slots written, the rest of the data as it is */
  posesOnly = false,
) {
  const tSegment = (values as { tSegment?: ValueCell<{ array: Float32Array }> }).tSegment?.ref.value;
  const layout = CulledData.layout(
    getTransformParts(values, n),
    values.aInstance.ref.value as Uint32Array,
    segments,
    getCulledSegmentData(segments, tSegment),
    n,
    MergedCullerLimits.dataBlockInstances,
  );
  ensure(data, layout.vec4s * 16);
  const buffer = data.ref.value;
  if (!posesOnly) {
    layout.write((d, o) => buffer.write(d, o), { poses: false, runs: false });
    runWriter.write(buffer, layout.runs, dispatchItems);
  }
  expander.expand(layout, buffer, dispatchItems);
}

const EmptyInt = new Int32Array(0);
const EmptyUint = new Uint32Array(0);

export function createMergedCuller(
  compute: ComputeContext,
  /** for a texture bound while there is no pyramid */
  resources: GraphicsResources,
  merged: Merged,
  members: readonly MergeableValues[],
  indexed: boolean,
  /**
   * items drawn by the shared indices of their pattern (no elements): an indexed record's first index 0, its base
   * vertex its segment's first item's (`ItemIndexPattern`), so the indices cover the largest segment's items alone
   */
  items?: ItemIndexPattern,
): MergedCuller {
  const { values, segments } = merged;
  /**
   * the members' cells the layout reads, resolved once into objects of one shape: a member's values object's named
   * reads are slow (thousands of members a layout, `getSegmentData`, `layoutRegions`, `layoutCells`)
   */
  const layoutCellsOf = members.map((m) => ({
    drawCount: m.drawCount as ValueCell<number>,
    lodLevels: m.lodLevels as ValueCell<LodLevelsValue> | undefined,
    invariantBoundingSphere: m.invariantBoundingSphere as ValueCell<Sphere3D>,
    boundingSphere: m.boundingSphere as ValueCell<Sphere3D>,
    padding: (m as { padding?: ValueCell<number> }).padding,
    instanceCellSize: m.instanceCellSize as ValueCell<number>,
  }));
  /** segment `s`'s first slot in a region (`MergedSegments.regionBases`, its data's first without) */
  const regionBaseOf = (s: number) => (segments.regionBases ?? segments.instanceBases)[s];
  const stride = indexed ? IndirectRecordBytes.indexed : IndirectRecordBytes.vertex;
  const words = stride / 4;

  const createBuffer = (bytes: number, label: string) => compute.createBuffer(Math.max(16, bytes), label);
  const buffer = (label: string) => ValueCell.create(createBuffer(16, label));
  const instanceCells = buffer('culled instance cells');
  const segmentData = buffer('culled segment data');
  const cellSpheres = buffer('culled cell spheres');
  const cellLevels = buffer('culled cell levels');
  const counters = buffer('culled counters');
  const args = buffer('culled draw arguments');
  const secondArgs = buffer('culled draw arguments, second phase');
  /** the regions' instance indices, then the identity region's (`identityBase`) */
  const indices = buffer('culled instance indices');
  /** the instances' data (`CulledData`) the indices index */
  const data = buffer('culled instance data');
  const culledValues: CulledValues = {
    sCulledIndices: indices,
    sCulledData: data,
    dCulled: ValueCell.create(true),
    dCulledPaged: ValueCell.create(false),
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
    ...createHiZValues(createEmptyHiZ(resources)),
    uClipObjectCount: ValueCell.create(0),
    uClipBase: ValueCell.create(0),
    uObjectTransform: ValueCell.create(Mat4.identity()),
    uObjectScale: ValueCell.create(1),
    uGhostRegion: ValueCell.create(0),
    segmentData,
    cellSpheres,
    cellLevels,
  };
  const emptyHiZ = cellValues.tHiZ.ref.value;
  // the sphere tests' uniforms shared with the cells' kernel (the same cells, set once per dispatch)
  const { uCellCount: _, cellSpheres: __, cellLevels: ___, ...sphereTestValues } = cellValues;
  const copyValues = {
    ...sphereTestValues,
    uInstanceCount: ValueCell.create(0),
    uSegmentCount: ValueCell.create(0),
    uRasterRadius: ValueCell.create(0),
    uRasterScale: ValueCell.create(0),
    uRasterNear: ValueCell.create(0),
    uRasterLevels: ValueCell.create(0),
    uClipInstance: ValueCell.create(0),
    uPhase: ValueCell.create(0),
    instanceCells,
    cellLevels,
    counters,
    indices,
    data,
  };
  const argsValues = {
    uRecordCount: ValueCell.create(0),
    uRecordWords: ValueCell.create(words),
    uPhase: copyValues.uPhase,
    counters,
    args,
    second: secondArgs,
  };
  let kernels: ComputeKernel[] | undefined;
  let cellBuilder: SceneCellBuilder | undefined;
  /** created at the first instances' data (its kernel at once: the data complete before the first draw) */
  let poseExpander: CulledPoseExpander | undefined;
  /** created with the first instances' data or identity region, as `poseExpander` */
  let runWriter: U32RunWriter | undefined;
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
    cellBuilder = createSceneCellBuilder(compute);
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
  let chunkCapacity = 0;
  /** whether the layout has raster regions (from the first dispatch with raster parameters on) */
  let rasterRegions = false;
  /** whether the layout has the ghost's region, the last (from the first dispatch with the ghost on) */
  let ghostRegion = false;
  let ghostDispatched = false;
  /** the regions of the layout: the levels', the raster's, the ghost's */
  let regionCount = 1;
  /** the ghost's region, 0 for none */
  let ghostIndex = 0;
  /** the ghost's list (`getGhostLists`), for the levels' lists and the records it was built with */
  let ghostLists: MultiDrawBaseData[] | undefined;
  let ghostListsLevels: readonly MultiDrawBaseData[] | undefined;
  let ghostListsArgs: StorageBuffer | undefined;
  /** `getSegmentData` of the layout */
  let segmentDataCopy: Float32Array | undefined;
  /** the largest scale of each segment's instances (`rasterized` scales the largest sphere by the instance's) */
  let segmentScales = new Float32Array(0);
  let rasterDispatched = false;
  let secondPending = false;
  let secondDispatched = false;
  /**
   * What the cells were split from (`syncLayout`), per member: its grid, the grid's layout and placement
   * (`moveInstanceGrid`), its transforms' and bounds' versions. A member's grid moved from the same layout builds the
   * cells again (`moveLayout`), not the rest of the layout
   */
  const splitGrids: (InstanceGrid | undefined)[] = [];
  const splitLayouts: (InstanceGrid | undefined)[] = [];
  const splitPlacements: Mat4[] = [];
  const splitTransforms: number[] = [];
  const splitBounds: number[] = [];
  /** the cells' grids (`layoutCells`); the cells built on the GPU at the next dispatch where pending */
  let sceneGrid: SceneGrid | undefined;
  let cellsPending = false;
  /** `segmentScales` of the split's instances, what a move scales */
  let splitScales = new Float32Array(0);
  let scalesSplit = false;

  /**
   * Per segment what the layout reads of its member (`refreshSegments`): its segment data's row but the first slot
   * (`getSegmentData`), its draws and its levels' counts (`layoutRegions`), its bounds and invariant radius
   * (`layoutCells`)
   */
  let segmentRows = new Float32Array(0);
  let segmentDraws = new Uint32Array(0);
  /** per segment and level (`MaxCullLevels` a segment) its count, where `segmentLeveled` (else its draws at each) */
  let segmentLevelCounts = new Uint32Array(0);
  let segmentLeveled = new Uint8Array(0);
  const segmentSpheres: Sphere3D[] = [];
  let segmentRadii = new Float32Array(0);

  function refreshSegment(s: number) {
    const member = layoutCellsOf[s];
    const o = s * 12;
    const { center, radius } = member.invariantBoundingSphere.ref.value;
    segmentRows[o] = center[0];
    segmentRows[o + 1] = center[1];
    segmentRows[o + 2] = center[2];
    segmentRows[o + 3] = radius;
    segmentRows[o + 5] = member.padding?.ref.value ?? 0;
    // gridless members draw at full detail at each level they overlap, as the CPU culling (`cullSegment`)
    const gridless = !hasInstanceCells(member);
    const memberLevels = !gridless ? member.lodLevels?.ref.value : undefined;
    for (let j = 0; j < 4; ++j) segmentRows[o + 8 + j] = memberLevels?.[j]?.[3] ?? 1;
    segmentDraws[s] = member.drawCount.ref.value;
    const leveled = levelCount > 0 && memberLevels !== undefined;
    segmentLeveled[s] = leveled ? 1 : 0;
    if (leveled) {
      for (let j = 0; j < levelCount; ++j) segmentLevelCounts[s * MaxCullLevels + j] = memberLevels[j][2];
    }
    segmentSpheres[s] = member.boundingSphere.ref.value;
    segmentRadii[s] = radius;
  }

  /** the segments' reads (`refreshSegment`) */
  function refreshSegments() {
    if (segmentDraws.length !== segmentCount) {
      segmentRows = new Float32Array(segmentCount * 12);
      segmentDraws = new Uint32Array(segmentCount);
      segmentLevelCounts = new Uint32Array(segmentCount * MaxCullLevels);
      segmentLeveled = new Uint8Array(segmentCount);
      segmentRadii = new Float32Array(segmentCount);
      segmentSpheres.length = segmentCount;
    }
    for (let s = 0; s < segmentCount; ++s) refreshSegment(s);
  }

  function syncLayout() {
    MergedCullerCounts.layouts += 1;
    scalesSplit = true;
    const lodLevels = members[0].lodLevels?.ref.value as LodLevelsValue | undefined;
    levelCount = lodLevels ? Math.min(lodLevels.length, MaxCullLevels) : 0;
    instanceCount = values.instanceCount.ref.value;
    segmentCount = segments.count;
    refreshSegments();

    setSegmentScales();
    // what a move or poses alone are told from (`isLayoutMoved`, `isPosesChange`)
    for (let s = 0; s < segmentCount; ++s) {
      const member = members[s];
      const grid = member.instanceGrid.ref.value;
      splitGrids[s] = grid;
      splitLayouts[s] = grid.layout ?? grid;
      splitPlacements[s] = Mat4.copy(splitPlacements[s] ?? Mat4(), grid.placement ?? IdentityPlacement);
      splitTransforms[s] = member.aTransform.ref.version;
      splitBounds[s] = member.boundingSphere.ref.version;
    }
    splitGrids.length = splitLayouts.length = splitPlacements.length = segmentCount;
    splitTransforms.length = splitBounds.length = segmentCount;

    // the clip objects after the segments (`syncClip`)
    ensure(segmentData, segmentCount * 48 + Clip.MaxObjectCount * ClipStorageVec4s * 16);
    segmentDataCopy = getSegmentData();
    segmentData.ref.value.write(segmentDataCopy);
    syncClip(true);
    layoutCells();

    ValueCell.updateIfChanged(levelCountCell, levelCount);
    ValueCell.updateIfChanged(copyValues.uInstanceCount, instanceCount);
    ValueCell.updateIfChanged(copyValues.uSegmentCount, segmentCount);
    layoutRegions();
  }

  /**
   * The cells: one grid over all instances, by position (`getSceneGrid`, from the members' bounds and spheres by their
   * segments' largest scales), the cells built from the instances' data on the GPU at the next dispatch
   * (`SceneCellBuilder`: the data written by then, `syncInstances`)
   */
  function layoutCells() {
    const radii = new Float32Array(segmentCount);
    for (let s = 0; s < segmentCount; ++s) {
      // read again: a move changes their bounds alone
      const member = layoutCellsOf[s];
      segmentSpheres[s] = member.boundingSphere.ref.value;
      segmentRadii[s] = member.invariantBoundingSphere.ref.value.radius;
      radii[s] = segmentRadii[s] * segmentScales[s];
    }
    sceneGrid = getSceneGrid(
      segmentSpheres,
      radii,
      segments.instanceCounts,
      instanceCount,
      MergedCullerLimits.sceneCellInstances,
    );
    cellCount = sceneGrid.cellCount;
    ensure(instanceCells, instanceCount * 8);
    ensure(cellSpheres, cellCount * 16);
    ensure(cellLevels, cellCount * 4);
    ValueCell.updateIfChanged(cellValues.uCellCount, cellCount);
    cellsPending = true;
  }

  /**
   * The regions of the layout and their records, apart from the cells (`syncLayout`): the raster's and the ghost's
   * turned on at a dispatch lay these out again, not the cells nor the instances' data
   */
  function layoutRegions() {
    const levels = Math.max(1, levelCount);
    // the levels' regions, with raster regions those again, then the ghost's
    const regions = (rasterRegions ? 2 : 1) * levels + (ghostRegion ? 1 : 0);
    regionCount = regions;
    ghostIndex = ghostRegion ? regions - 1 : 0;
    identityBase = regions * instanceCount;

    ensure(counters, regions * segmentCount * 4);
    ensure(args, regions * segmentCount * stride);
    ensure(secondArgs, regions * segmentCount * stride);
    ensure(indices, (regions + 1) * instanceCount * 4);
    // the identity region: every instance, in its own order or the one set (`setIdentityOrder`)
    writeIdentityRegion();

    // the records' words but the instance count: the level's count, first vertex or index, the region's first instance
    const records = new Uint32Array(regions * segmentCount * words);
    chunkCapacity = 0;
    // per segment its draws and levels as read (`refreshSegments`)
    for (let k = 0; k < regions; ++k) {
      // a raster region as its level's, the ghost's as the coarsest
      const j = ghostRegion && k === ghostIndex ? levels - 1 : k % levels;
      for (let s = 0; s < segmentCount; ++s) {
        const o = (k * segmentCount + s) * words;
        records[o] = segmentLeveled[s] ? segmentLevelCounts[s * MaxCullLevels + j] : segmentDraws[s];
        // quads of 6 vertices (or indices) a sphere
        if (k === j && !(ghostRegion && k === ghostIndex))
          chunkCapacity += Math.ceil(records[o] / 6 / 64) * segments.instanceCounts[s];
        const first = k * instanceCount + regionBaseOf(s);
        if (indexed && items) {
          // the first index (`elementOffsets` in bytes) as the segment's first item, by its first vertex
          records[o + 3] = (segments.elementOffsets[s] / 4 / items.pattern.length) * items.verticesPerItem;
          records[o + 4] = first;
        } else if (indexed) {
          records[o + 2] = segments.elementOffsets[s] / 4;
          records[o + 4] = first;
        } else {
          records[o + 2] = segments.drawBases[s];
          records[o + 3] = first;
        }
      }
    }
    args.ref.value.write(records);
    secondArgs.ref.value.write(records);

    ValueCell.updateIfChanged(argsValues.uRecordCount, regions * segmentCount);
    ghostLists = undefined;
    ValueCell.updateIfChanged(copyValues.uRasterLevels, rasterRegions ? levels : 0);
  }

  /** `setIdentityOrder`'s order, for as many instances as the layout's */
  let identityOrder: Uint32Array | undefined;
  /** the identity region's run written by a dispatch not submitted yet (a CPU write after it submits it first) */
  let identityPending = false;
  function writeIdentityRegion() {
    if (identityOrder && identityOrder.length !== instanceCount) identityOrder = undefined;
    if (identityOrder) {
      if (identityPending) compute.submit();
      identityPending = false;
      indices.ref.value.write(identityOrder, identityBase * 4);
      return;
    }
    // the instances' own order a run, written on the GPU (written again where the regions grow: no upload)
    runWriter ??= createU32RunWriter(compute);
    runWriter.write(
      indices.ref.value,
      [{ word: identityBase, count: instanceCount, start: 0, step: 1, bits: 0 }],
      dispatchItems,
      false,
    );
    identityPending = true;
  }

  /**
   * Whether the layout's split holds for the members as they are: each grid of the layout it was split from (the same
   * grid with the same transforms, or one moved from it, `moveInstanceGrid`), all else of the layout key the same; the
   * members moved alike or not (the cells built again)
   */
  function isLayoutMoved() {
    if (!isSameKey(getFixedKey(), fixedKey) || segments.count !== splitGrids.length) return false;
    for (let s = 0; s < segmentCount; ++s) {
      const member = members[s];
      const grid = member.instanceGrid.ref.value;
      if ((grid.layout ?? grid) !== splitLayouts[s]) return false;
      // the grid as split, its transforms or bounds (a gridless member's cell) changed: not a move
      if (
        grid === splitGrids[s] &&
        (member.aTransform.ref.version !== splitTransforms[s] || member.boundingSphere.ref.version !== splitBounds[s])
      ) {
        return false;
      }
    }
    return true;
  }
  const moveSegment = Mat4();

  /**
   * Whether the members' poses changed alone since the layout (`Particle.withPoses`'s
   * frames): each grid the one split and pending (its instances in their own order, the slots' words theirs), all
   * else of the layout key the same; a member with a built grid changed is laid out again
   */
  function isPosesChange() {
    if (!isSameKey(getFixedKey(), fixedKey) || segments.count !== splitGrids.length) return false;
    for (let s = 0; s < segmentCount; ++s) {
      const member = members[s];
      const grid = member.instanceGrid.ref.value;
      if (grid !== splitGrids[s]) return false;
      if (!grid.pending && member.aTransform.ref.version !== splitTransforms[s]) return false;
    }
    return true;
  }

  /** The poses changed alone (`isPosesChange`): the cells built again from the instances as they are now */
  function posesLayout() {
    MergedCullerCounts.poses += 1;
    for (let s = 0; s < segmentCount; ++s) {
      splitTransforms[s] = members[s].aTransform.ref.version;
      splitBounds[s] = members[s].boundingSphere.ref.version;
    }
    layoutCells();
  }
  /** segment `s`'s move since the cells were laid out (`placement · placement_split⁻¹`) into `out` */
  function getSegmentMove(s: number, out: Mat4) {
    const placement = members[s].instanceGrid.ref.value.placement ?? IdentityPlacement;
    return Mat4.mul(out, placement, Mat4.invert(moveScratch, splitPlacements[s]));
  }

  /**
   * The segments moved since the layout (`isLayoutMoved`): their scales by their moves' largest scales, the cells built
   * again from the moved instances (`layoutCells`); not laid out again, the segments' data and the records kept (a
   * placement moved as a whole, the release of a drag)
   */
  function moveLayout() {
    MergedCullerCounts.moves += 1;
    for (let s = 0; s < segmentCount; ++s) {
      segmentScales[s] = splitScales[s] * Mat4.getMaxScaleOnAxis(getSegmentMove(s, moveSegment));
    }
    layoutCells();
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

  /**
   * per segment its invariant sphere, its first instance and largest sphere's radius (`padding`: the largest size times
   * the size factor), its level scales (as `tSegment`'s) (`segmentData`): the rows as read (`refreshSegments`), the
   * first instances as laid out
   */
  function getSegmentData() {
    const out = segmentRows.slice(0, segmentCount * 12);
    const bases = new Uint32Array(out.buffer);
    for (let s = 0; s < segmentCount; ++s) bases[s * 12 + 4] = regionBaseOf(s);
    return out;
  }

  /**
   * The clip objects cells and instances are culled by where they clip them whole (`getCullClipCount`), instances' centers
   * clipped too with instance clipping (`getInstanceClipCount`, the vertex stage's), after the segments in
   * `segmentData`; where they changed (or `force`, the buffer written anew) written, true
   */
  const clipVersions = [-1, -1, -1, -1];
  function syncClip(force: boolean) {
    const count = getCullClipCount(values);
    const frame = values.uClipObjectFrame.ref.version;
    const shape = values.uClipObjectShape.ref.version;
    const instance = getInstanceClipCount(values) > 0 ? 1 : 0;
    const v = clipVersions;
    if (!force && v[0] === count && v[1] === frame && v[2] === shape && v[3] === instance) return false;
    v[0] = count;
    v[1] = frame;
    v[2] = shape;
    v[3] = instance;
    if (count > 0) {
      const data = new Float32Array(count * ClipStorageVec4s * 4);
      packStorageClipObjects(values, data);
      segmentData.ref.value.write(data, segmentCount * 48);
    }
    ValueCell.updateIfChanged(cellValues.uClipObjectCount, count);
    ValueCell.updateIfChanged(cellValues.uClipBase, segmentCount * 3);
    ValueCell.updateIfChanged(copyValues.uClipInstance, instance);
    return true;
  }

  /**
   * the instances' data (`CulledData`), what the kernels and the culled programs read by index; `moved`: the segments'
   * scales were set by the move (`moveLayout`)
   */
  function syncInstances(moved = false, posesOnly = false) {
    if (!moved && !posesOnly) setSegmentScales();
    if (scalesSplit) {
      splitScales = segmentScales.slice();
      scalesSplit = false;
    }
    poseExpander ??= createCulledPoseExpander(compute);
    runWriter ??= createU32RunWriter(compute);
    writeCulledData(data, ensure, values, segments, instanceCount, poseExpander, runWriter, posesOnly);
  }

  /** debug mode: what the raster split classifies by (`rasterized`), logged where it changed */
  let rasterInputsLogged = '';
  function logRasterInputs(raster: MergedRasterParams) {
    const data = getSegmentData();
    const range = (o: number) => {
      let min = Infinity;
      let max = -Infinity;
      for (let s = 0; s < segmentCount; ++s) {
        min = Math.min(min, data[s * 12 + o]);
        max = Math.max(max, data[s * 12 + o]);
      }
      return `${min.toFixed(2)}-${max.toFixed(2)}`;
    };
    const levels: string[] = [];
    for (let j = 0; j < Math.max(1, levelCount); ++j) {
      levels.push(`${levelRanges[j * 4].toFixed(1)}-${levelRanges[j * 4 + 1].toFixed(1)} x${range(8 + j)}`);
    }
    const text =
      `sphere raster split: ${JSON.stringify(raster)}, ${segmentCount} segments, invariant radius ${range(3)}, ` +
      `largest sphere ${range(5)}, levels (distances, scales) ${levels.join('; ')}`;
    if (text !== rasterInputsLogged) console.log(text);
    rasterInputsLogged = text;
  }

  // the keys as numbers in arrays kept, compared in place: a sync each update of the object (a clip's drag: each frame)
  // over thousands of members (the yeast cell's 5,639), not strings built of them all (1.2 ms a sync)
  let fixedScratch: Float64Array = new Float64Array(0);
  let layoutScratch: Float64Array = new Float64Array(0);
  const keyArray = (array: Float64Array, length: number): Float64Array =>
    array.length === length ? array : new Float64Array(length);
  /**
   * the layout key but what a move changes (the grids, the transforms, the bounds): what a moved layout keeps; written
   * into `out` (another array where its length differs)
   */
  const writeFixedKey = (out: Float64Array) => {
    const key = keyArray(out, 4 + members.length * 3);
    key[0] = values.instanceCount.ref.value;
    key[1] = segments.count;
    key[2] = rasterRegions ? 1 : 0;
    key[3] = ghostRegion ? 1 : 0;
    for (let i = 0, o = 4, il = members.length; i < il; ++i, o += 3) {
      const m = members[i];
      key[o] = m.invariantBoundingSphere.ref.version;
      key[o + 1] = m.drawCount.ref.value;
      key[o + 2] = m.lodLevels?.ref.version ?? -1;
    }
    return key;
  };
  /** what the cells and records are written from: the grids and the instances they split, the counts and the levels */
  const writeLayoutKey = (out: Float64Array) => {
    const fixed = 4 + members.length * 3;
    const key = keyArray(out, fixed + members.length * 3);
    fixedScratch = writeFixedKey(fixedScratch);
    key.set(fixedScratch);
    for (let i = 0, o = fixed, il = members.length; i < il; ++i, o += 3) {
      const m = members[i];
      key[o] = m.instanceGrid.ref.version;
      key[o + 1] = m.aTransform.ref.version;
      key[o + 2] = m.boundingSphere.ref.version;
    }
    return key;
  };
  const isSameKey = (a: Float64Array, b: Float64Array) => {
    if (a.length !== b.length) return false;
    for (let i = 0, il = a.length; i < il; ++i) if (a[i] !== b[i]) return false;
    return true;
  };
  const getFixedKey = () => (fixedScratch = writeFixedKey(fixedScratch));
  const getLayoutKey = () => (layoutScratch = writeLayoutKey(layoutScratch));
  let fixedKey: Float64Array = new Float64Array(0);
  const getInstanceVersions = () =>
    `${values.aTransform.ref.version},${values.instancePoses?.ref.version},${values.aInstance.ref.version},${values.aSegment.ref.version}`;
  /** none before the first layout */
  let layoutKey: Float64Array = new Float64Array(0);
  let instanceVersions = '';
  let transformVersion = -1;
  /** advanced where the layout or the instances changed: a dispatch culled with an older one is not current */
  let dataVersion = 0;
  function sync() {
    const key = getLayoutKey();
    const layoutChanged = !isSameKey(key, layoutKey);
    const laidOut = layoutKey.length > 0;
    // a placement moved as a whole (grids moved from those split): the cells moved, not split again
    const moved = layoutChanged && laidOut && isLayoutMoved();
    // the poses changed alone: the cells built again, the regions and the rest of the data kept
    const posesOnly = layoutChanged && !moved && laidOut && isPosesChange();
    if (layoutChanged) {
      layoutKey = key.slice();
      if (moved) moveLayout();
      else if (posesOnly) posesLayout();
      else {
        fixedKey = getFixedKey().slice();
        syncLayout();
      }
    }
    // the object's transient transform (a drag's frame): two uniforms the kernels move each sphere by, not the layout
    const transform = values.uObjectTransform;
    if (transform.ref.version !== transformVersion) {
      transformVersion = transform.ref.version;
      ValueCell.update(
        cellValues.uObjectTransform,
        Mat4.copy(cellValues.uObjectTransform.ref.value, transform.ref.value),
      );
      ValueCell.update(cellValues.uObjectScale, Mat4.getMaxScaleOnAxis(transform.ref.value));
      dataVersion += 1;
    }
    const versions = getInstanceVersions();
    if (layoutChanged || versions !== instanceVersions) {
      instanceVersions = versions;
      syncInstances(moved, posesOnly);
      dataVersion += 1;
    }
  }
  /**
   * A region turned on (`turnOn`: the raster's, the ghost's): the regions laid out again (`layoutRegions`) where the
   * layout is otherwise current, not the cells nor the instances' data; else all of it (`sync`)
   */
  function setRegions(turnOn: () => void) {
    const current = isSameKey(getLayoutKey(), layoutKey);
    turnOn();
    if (!current) {
      sync();
      return;
    }
    layoutKey = getLayoutKey().slice();
    fixedKey = getFixedKey().slice();
    layoutRegions();
    dataVersion += 1;
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
   * pyramid's camera; while they are the same its regions' indices and records are current (a still camera culls
   * nothing again, once tested twice against a pyramid of its own view: the first such pyramid is of a frame culled
   * against another view's, which can lack occluders and occlude less; a round counts a new build only, a pyramid
   * served again while the Hi-Z pass skips its builds is still of the frame it was built in)
   */
  const dispatched = new Float64Array(7 * 4 + 2).fill(NaN);
  const dispatchedLevels = Mat4.zero().fill(NaN);
  let dispatchedVersion = -1;
  /** the pyramid's camera: its view and projection, its second view's (`OcclusionPyramid.second`), the split */
  const dispatchedPyramid = new Float64Array(65).fill(NaN);
  const setPyramidCamera = (out: Float64Array, pyramid: OcclusionPyramid) => {
    out.set(pyramid.view, 0);
    out.set(pyramid.projection, 16);
    const { second } = pyramid;
    if (second) {
      out.set(second.view, 32);
      out.set(second.projection, 48);
    } else out.fill(0, 32, 64);
    out[64] = second ? second.split : 1;
  };
  const pyramidCamera = new Float64Array(65);
  let dispatchedHiZ = false;
  let dispatchedPyramidVersion = -1;
  /** dispatches in a row against pyramids of the same camera, each a new build */
  let pyramidRounds = 0;
  function isPyramidCamera(pyramid: OcclusionPyramid) {
    setPyramidCamera(pyramidCamera, pyramid);
    for (let i = 0; i < 65; ++i) {
      if (dispatchedPyramid[i] !== pyramidCamera[i]) return false;
    }
    return true;
  }
  function isPyramidDispatched(pyramid: OcclusionPyramid | undefined) {
    if (!pyramid) return !dispatchedHiZ;
    return dispatchedHiZ && pyramidRounds >= 2 && isPyramidCamera(pyramid);
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

  /** the levels' lists, then the raster regions' where there are; `levelLists` the levels' of them (the same lists) */
  let lists: MultiDrawBaseData[] = [];
  let levelLists: MultiDrawBaseData[] = [];
  let listsLevels: readonly MultiDrawBaseData[] | undefined;
  let listsKey = '';

  /**
   * A still camera culls nothing again (`isDispatched`), so its records stay: read back once, each list then draws only
   * those with instances. A draw costs whether it draws or not, outside the GPU's timers: with thousands of units and
   * levels, most of them empty (the raster's instances are in its regions), the larger part of a frame. Any dispatch
   * or new lists restore all records
   */
  let recordsGeneration = 0;
  let recordsRead = -1;
  let recordsCompact = false;
  let recordsVersion = 0;
  let dispatchVersion = 0;
  /**
   * The records a dispatch can fill, per list (the levels', then the raster regions') and segment: without those the
   * segment's bounds rule out, so that a moving camera does not issue them either (`predictRecords`)
   */
  let predicted = new Uint8Array(0);
  /** each segment's bounds with the object's transient transform (`getDrawnBoundingSphere`, cached per segment) */
  const segmentBounds: DrawnBoundsValues[] = members.map((m) => ({
    boundingSphere: m.boundingSphere,
    uObjectTransform: values.uObjectTransform,
  }));
  let predictedNext = new Uint8Array(0);
  function setListRecords() {
    recordsVersion += 1;
    for (let k = 0; k < lists.length; ++k) {
      const l = lists[k];
      const records = new Uint32Array(segmentCount);
      let n = 0;
      for (let s = 0; s < segmentCount; ++s) if (predicted[k * segmentCount + s]) records[n++] = s;
      l.count = n;
      l.indirect = { ...l.indirect!, records: records.subarray(0, n), recordsVersion };
    }
  }
  function setSegmentScales() {
    if (segmentScales.length !== segmentCount) segmentScales = new Float32Array(segmentCount);
    for (let s = 0; s < segmentCount; ++s) {
      // the member's own transforms (the merged object's are its members' as parts), by their encoding
      segmentScales[s] = InstanceTransforms.getMaxScaleOfAll(
        getObjectTransforms(members[s] as unknown as TransformValues),
      );
    }
  }
  /**
   * From each segment's bounds, which hold its instances' spheres: conservative (a margin for the kernels' 32-bit
   * floats), the kernels decide per cell and instance. A record is left out where the bounds lie outside the frustum or
   * the level's range, and a level's where all the instances go to the raster (`rasterized` holds at their nearest
   * point, and so further away). False where nothing changed
   */
  function predictRecords(cameraPlane: Plane3D, frustum: Frustum3D, raster: MergedRasterParams | undefined) {
    const levels = Math.max(1, levelCount);
    const size = (rasterRegions ? 2 : 1) * levels * segmentCount;
    if (predictedNext.length !== size) predictedNext = new Uint8Array(size);
    const next = predictedNext.fill(0);
    const segmentInfo = segmentDataCopy!;
    const objectScale = cellValues.uObjectScale.ref.value;
    const c = [0, 0, 0];
    for (let s = 0; s < segmentCount; ++s) {
      // as drawn: moved by the object's transient transform
      let radius: number;
      {
        const drawn = getDrawnBoundingSphere(segmentBounds[s]);
        c[0] = drawn.center[0];
        c[1] = drawn.center[1];
        c[2] = drawn.center[2];
        radius = drawn.radius;
      }
      const r = radius * 1.001 + 1e-3;
      let inside = true;
      for (let i = 0; i < 6 && inside; ++i) {
        const { normal, constant } = frustum[i as 0];
        inside = normal[0] * c[0] + normal[1] * c[1] + normal[2] * c[2] + constant >= -r;
      }
      if (!inside) continue;
      const n = cameraPlane.normal;
      const d = n[0] * c[0] + n[1] * c[1] + n[2] * c[2] + cameraPlane.constant;
      for (let j = 0; j < levels; ++j) {
        const min = levelCount ? levelRanges[j * 4] : -Infinity;
        const max = levelCount ? levelRanges[j * 4 + 1] : Infinity;
        if (d + r < min || d - r > max) continue;
        let toRaster = false;
        if (raster) {
          next[(levels + j) * segmentCount + s] = 1;
          const o = s * 12;
          const largest =
            segmentInfo[o + 5] * segmentScales[s] * objectScale * segmentInfo[o + 8 + Math.min(j, 3)] * 1.001;
          const nearest = Math.max(d - r, min);
          toRaster =
            segmentInfo[o + 3] > 0 &&
            nearest - largest > raster.near &&
            largest * raster.scale <= raster.radius * Math.sqrt(Math.max(nearest * nearest - largest * largest, 0));
        }
        if (!toRaster) next[j * segmentCount + s] = 1;
      }
    }
    let same = next.length === predicted.length;
    for (let i = 0; same && i < next.length; ++i) same = next[i] === predicted[i];
    if (same) return false;
    predictedNext = predicted;
    predicted = next;
    return true;
  }
  function restoreRecords() {
    recordsGeneration += 1;
    if (!recordsCompact) return;
    recordsCompact = false;
    setListRecords();
  }
  function compactRecords() {
    if (recordsRead === recordsGeneration || recordsCompact || lists.length === 0) return;
    recordsRead = recordsGeneration;
    const generation = recordsGeneration;
    args.ref.value
      .readAsync()
      .then((data) => {
        if (generation !== recordsGeneration) return;
        const words = new Uint32Array(data);
        for (const l of lists) {
          const records = new Uint32Array(segmentCount);
          let n = 0;
          const base = l.indirect!.byteOffset / 4;
          for (let s = 0; s < segmentCount; ++s) {
            if (words[base + (s * stride) / 4 + 1] > 0) records[n++] = s;
          }
          l.count = n;
          l.indirect = { ...l.indirect!, records: records.subarray(0, n), recordsVersion: recordsVersion + 1 };
        }
        recordsCompact = true;
        recordsVersion += 1;
      })
      // the full lists stay: the same image
      .catch(() => {});
  }

  /** the uniforms of the ghost's list: the coarsest level at every distance */
  const ghostLod = ValueCell.create(Vec4());
  const ghostLodLevel = ValueCell.create(0);
  function getGhostList(levels: readonly MultiDrawBaseData[]): MultiDrawBaseData {
    const j = Math.max(1, levelCount) - 1;
    const level = levels[Math.min(j, levels.length - 1)];
    const lodLevels = members[0].lodLevels?.ref.value as LodLevelsValue | undefined;
    let uniforms: MultiDrawBaseData['uniforms'] = [];
    if (levelCount > 0 && lodLevels && level?.uniforms.length) {
      ValueCell.update(ghostLod, Vec4.set(ghostLod.ref.value, 0, Number.MAX_VALUE, 0, lodLevels[j][3]));
      ValueCell.updateIfChanged(ghostLodLevel, Math.min(j, 3));
      uniforms = level.uniforms.map(([name, cell]) =>
        name === 'uLod' ? ['uLod', ghostLod] : name === 'uLodLevel' ? ['uLodLevel', ghostLodLevel] : [name, cell],
      );
    }
    return {
      firsts: EmptyInt,
      counts: EmptyInt,
      offsets: EmptyInt,
      instanceCounts: EmptyInt,
      baseVertices: EmptyInt,
      baseInstances: EmptyUint,
      count: segmentCount,
      uniforms,
      indirect: { buffer: args.ref.value, byteOffset: ghostIndex * segmentCount * stride, indexed },
    };
  }

  function buildLists(levels: readonly MultiDrawBaseData[], key: string) {
    restoreRecords();
    listsLevels = levels;
    listsKey = key;
    const n = Math.max(1, levelCount);
    const list = (level: MultiDrawBaseData, k: number): MultiDrawBaseData => ({
      firsts: EmptyInt,
      counts: EmptyInt,
      offsets: EmptyInt,
      instanceCounts: EmptyInt,
      baseVertices: EmptyInt,
      baseInstances: EmptyUint,
      count: segmentCount,
      uniforms: level.uniforms,
      indirect: { buffer: args.ref.value, byteOffset: k * segmentCount * stride, indexed },
    });
    levelLists = levels.slice(0, n).map((level, j) => list(level, j));
    lists = rasterRegions ? [...levelLists, ...levelLists.map((level, j) => list(level, n + j))] : levelLists;
    if (predicted.length === lists.length * segmentCount) setListRecords();
  }

  return {
    values: culledValues,
    get identityBase() {
      return identityBase;
    },
    setIdentityOrder: (order) => {
      identityOrder = order;
      if (instanceCount > 0) writeIdentityRegion();
    },
    sync,
    get dataVersion() {
      return dataVersion;
    },
    areKernelsReady: () => getKernels(),
    dispatch: (cameraPlane, frustum, lod, pyramid, raster, batch, ghost = false) => {
      if (instanceCount === 0 || cellCount === 0) return false;
      if (raster && !rasterRegions) setRegions(() => (rasterRegions = true));
      if (ghost && !ghostRegion) setRegions(() => (ghostRegion = true));
      const ghostChanged = ghost !== ghostDispatched;
      // clipped instances culled: a change of the clip objects culls again
      if (syncClip(false)) dataVersion += 1;
      secondPending = false;
      secondDispatched = false;
      const [min, max] = values.uLod.ref.value as Vec4;
      setLevelRanges(lod);
      const rasterChanged =
        copyValues.uRasterRadius.ref.value !== (raster?.radius ?? 0) ||
        copyValues.uRasterScale.ref.value !== (raster?.scale ?? 0) ||
        copyValues.uRasterNear.ref.value !== (raster?.near ?? 0);
      const viewChanged = !isDispatched(cameraPlane, frustum, min, max);
      if (!rasterChanged && !ghostChanged && !viewChanged && isPyramidDispatched(pyramid)) {
        compactRecords();
        return true;
      }
      if (!getKernels()) return false;
      restoreRecords();
      if (predictRecords(cameraPlane, frustum, raster)) setListRecords();
      dispatchVersion += 1;
      ValueCell.updateIfChanged(copyValues.uRasterRadius, raster?.radius ?? 0);
      ValueCell.updateIfChanged(copyValues.uRasterScale, raster?.scale ?? 0);
      ValueCell.updateIfChanged(copyValues.uRasterNear, raster?.near ?? 0);
      rasterDispatched = !!raster;
      ghostDispatched = ghost;
      ValueCell.updateIfChanged(cellValues.uGhostRegion, ghost ? ghostIndex : 0);
      if (isDebugMode && raster) logRasterInputs(raster);
      if (!pyramid) pyramidRounds = 0;
      // a view of its own: levels of an earlier one (a draw compiling, none of this view yet) left out what stays out
      // until culled again
      else if (viewChanged || !dispatchedHiZ || !isPyramidCamera(pyramid)) pyramidRounds = 1;
      else if (pyramid.version !== dispatchedPyramidVersion) pyramidRounds += 1;
      dispatchedHiZ = !!pyramid;
      if (pyramid) {
        dispatchedPyramidVersion = pyramid.version;
        setPyramidCamera(dispatchedPyramid, pyramid);
      }
      setHiZValues(cellValues, pyramid, emptyHiZ);
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

      const records = regionCount * segmentCount;
      const k = kernels!;
      // the cells built first where laid out since (`layoutCells`)
      const builder = cellsPending ? cellBuilder! : undefined;
      if (builder) {
        builder.set(sceneGrid!, { instanceCount, segmentData, data, instanceCells, cellSpheres });
        cellsPending = false;
      }
      const clear = () => {
        builder?.clear();
        counters.ref.value.clear(0, records * 4);
      };
      const dispatches = () => {
        builder?.dispatch(dispatchItems);
        dispatchItems(k[0], cellCount);
        dispatchItems(k[1], instanceCount);
        dispatchItems(k[2], records);
      };
      if (batch) {
        batch.clears.push(clear);
        batch.dispatches.push(dispatches);
      } else {
        clear();
        compute.beginComputePass({ label: 'merged cell culling' });
        dispatches();
        compute.endComputePass();
      }
      setDispatched(cameraPlane, frustum, min, max);
      // without a pyramid the first phase left out nothing occluded
      secondPending = !!pyramid;
      return true;
    },
    get secondPending() {
      return secondPending;
    },
    get ghostDispatched() {
      return ghostDispatched;
    },
    getGhostLists: (levels) => {
      if (!ghostLists || ghostListsLevels !== levels || ghostListsArgs !== args.ref.value) {
        ghostListsLevels = levels;
        ghostListsArgs = args.ref.value;
        ghostLists = [getGhostList(levels)];
      }
      return ghostLists;
    },
    dispatchSecond: (pyramid, batch) => {
      if (!secondPending) return false;
      secondPending = false;
      if (!pyramid || !kernels) return false;
      setHiZValues(cellValues, pyramid, emptyHiZ);
      const records = regionCount * segmentCount;
      const k = kernels;
      // the phase a uniform the kernels write when dispatched
      const dispatches = () => {
        ValueCell.updateIfChanged(copyValues.uPhase, 1);
        dispatchItems(k[0], cellCount);
        dispatchItems(k[1], instanceCount);
        dispatchItems(k[2], records);
        ValueCell.updateIfChanged(copyValues.uPhase, 0);
      };
      if (batch) {
        batch.dispatches.push(dispatches);
      } else {
        compute.beginComputePass({ label: 'merged cell culling, second phase' });
        dispatches();
        compute.endComputePass();
      }
      secondDispatched = true;
      return true;
    },
    get secondDispatched() {
      return secondDispatched;
    },
    getSecondLists: (levels, withRaster = true) => {
      const key = `${segmentCount},${levelCount},${args.ref.value.id},${rasterRegions}`;
      if (levels !== listsLevels || key !== listsKey) buildLists(levels, key);
      const buffer = secondArgs.ref.value;
      return (withRaster ? lists : levelLists).map((l) => ({ ...l, indirect: { ...l.indirect!, buffer } }));
    },
    getLists: (levels, withRaster = true) => {
      // the raster regions too where they exist and the raster does not draw them (another view, not ready); the
      // levels' lists the same objects either way, so that their draws (and recordings of them) stay the same
      const key = `${segmentCount},${levelCount},${args.ref.value.id},${rasterRegions}`;
      if (levels !== listsLevels || key !== listsKey) buildLists(levels, key);
      return withRaster ? lists : levelLists;
    },
    getRasterInput: () => {
      if (!rasterRegions || !rasterDispatched) return undefined;
      const levels = Math.max(1, levelCount);
      return {
        indices: indices.ref.value,
        data: data.ref.value,
        records: args.ref.value,
        secondRecords: secondArgs.ref.value,
        get secondDispatched() {
          return secondDispatched;
        },
        recordWords: words,
        levels,
        levelCount,
        lodLevels: members[0].lodLevels?.ref.value as LodLevelsValue | undefined,
        instanceCount,
        segmentCount,
        chunkCapacity,
        dispatchVersion,
        recordBase: levels * segmentCount,
        slotBase: levels * instanceCount,
      };
    },
    get recordsCompact() {
      return recordsCompact;
    },
    get recordsVersion() {
      return recordsVersion;
    },
    destroy: () => {
      if (kernels) for (const k of kernels) k.destroy();
      cellBuilder?.destroy();
      poseExpander?.destroy();
      runWriter?.destroy();
      emptyHiZ.destroy();
      const cells = [instanceCells, segmentData, cellSpheres, cellLevels, counters, args, secondArgs, indices, data];
      for (const cell of cells) cell.ref.value.destroy();
    },
  };
}
