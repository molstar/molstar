/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { OcclusionPyramid } from './compute/cell-cull.js';
import type { StorageBuffer } from './compute.js';
import type { RenderableValues, Values, RenderableSchema, BaseValues } from './renderable/schema.js';
import type {
  ComputeRenderItem,
  GlobalBindings,
  GraphicsProgram,
  GraphicsRenderItem,
  GraphicsRenderVariant,
  MultiDrawBaseData,
  Transparency,
} from './render-item.js';
import { ValueCell } from '@molstar/core/util/value-cell';
import { idFactory } from '@molstar/core/util/id-factory';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import type { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import type { GraphicsStats } from './stats.js';
import { isTimingMode } from '@molstar/core/util/debug';
import type { InstanceGrid } from '@molstar/core/math/geometry/instance-grid';
import { Clip } from '@molstar/core/util/clip';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const p3distanceToPoint = Plane3D.distanceToPoint;
const f3intersectsSphere3D = Frustum3D.intersectsSphere3D;
const s3fromArray = Sphere3D.fromArray;

const getNextRenderableId = idFactory();

// Opaque, unforgeable token identifying "a frame" for cull-cache reuse. Deliberately an
// object identity rather than a plain number: several independent owners (Canvas3D's main
// render loop, ImagePass, PickHelper, RayHelper, ...) each mint their own frame tokens while
// culling the SAME underlying renderables - if `frame` were a plain counter, two unrelated
// owners could coincidentally reach the same integer at the same time and wrongly reuse each
// other's cull result. Object identity (`===`) makes that structurally impossible: a token
// from one owner can never equal a token from another, no coordination between owners needed.
declare const frameBrand: unique symbol;
export type Frame = { readonly [frameBrand]: never };
/** Mint a new frame token - call this whenever an owner starts a new logical "frame". */
export function createFrame(): Frame {
  return {} as Frame;
}

export type RenderableState = {
  disposed: boolean;
  visible: boolean;
  pickable: boolean;
  colorOnly: boolean;
  opaque: boolean;
  writeDepth: boolean;
};

/** A state whose writes notify its observers (`RenderableState.observe`), as `ValueCell.observe` does for values */
class ObservableRenderableState implements RenderableState {
  observers: Set<ValueCell.Observer> | undefined = undefined;
  private _disposed: boolean;
  private _visible: boolean;
  private _pickable: boolean;
  private _colorOnly: boolean;
  private _opaque: boolean;
  private _writeDepth: boolean;

  constructor(state: RenderableState) {
    this._disposed = state.disposed;
    this._visible = state.visible;
    this._pickable = state.pickable;
    this._colorOnly = state.colorOnly;
    this._opaque = state.opaque;
    this._writeDepth = state.writeDepth;
  }

  private notify() {
    const { observers } = this;
    if (observers) for (const o of observers) o.changed = true;
  }

  get disposed() {
    return this._disposed;
  }
  set disposed(v: boolean) {
    if (v === this._disposed) return;
    this._disposed = v;
    this.notify();
  }
  get visible() {
    return this._visible;
  }
  set visible(v: boolean) {
    if (v === this._visible) return;
    this._visible = v;
    this.notify();
  }
  get pickable() {
    return this._pickable;
  }
  set pickable(v: boolean) {
    if (v === this._pickable) return;
    this._pickable = v;
    this.notify();
  }
  get colorOnly() {
    return this._colorOnly;
  }
  set colorOnly(v: boolean) {
    if (v === this._colorOnly) return;
    this._colorOnly = v;
    this.notify();
  }
  get opaque() {
    return this._opaque;
  }
  set opaque(v: boolean) {
    if (v === this._opaque) return;
    this._opaque = v;
    this.notify();
  }
  get writeDepth() {
    return this._writeDepth;
  }
  set writeDepth(v: boolean) {
    if (v === this._writeDepth) return;
    this._writeDepth = v;
    this.notify();
  }
}

export namespace RenderableState {
  /** A state with the values of `state` whose writes notify observers (`observe`) */
  export function create(state: RenderableState): RenderableState {
    return new ObservableRenderableState(state);
  }

  export function isObservable(state: RenderableState) {
    return state instanceof ObservableRenderableState;
  }

  /**
   * `observer.changed` is set when a field of `state` changes, until `unobserve`; states that are not observable
   * (`isObservable`) throw
   */
  export function observe(state: RenderableState, observer: ValueCell.Observer) {
    if (!(state instanceof ObservableRenderableState)) throw new Error('state not observable');
    (state.observers ??= new Set()).add(observer);
  }

  export function unobserve(state: RenderableState, observer: ValueCell.Observer) {
    if (state instanceof ObservableRenderableState) state.observers?.delete(observer);
  }
}

export interface Renderable<T extends RenderableValues> {
  readonly id: number;
  readonly materialId: number;
  /** For pipeline requests (`resources.preparePipelines`) */
  readonly renderItem: GraphicsRenderItem;
  readonly values: T;
  readonly state: RenderableState;

  cull: (
    cameraPlane: Plane3D,
    frustum: Frustum3D,
    isOccluded: ((s: Sphere3D) => boolean) | null,
    stats: GraphicsStats,
    frame: Frame,
    lod: LodView,
  ) => void;
  uncull: (lod: LodView) => void;
  cullSimple: (d: number, radius: number, scale: number, lod: LodView) => void;
  /**
   * Culls on the GPU, where the renderable does (merged render objects, see `merged-culler.ts`): called by
   * `Renderer.prepareCull` before the passes of `frame` (no pass open), its `cull` with the same planes then draws the
   * result
   */
  prepareCull?: (
    cameraPlane: Plane3D,
    frustum: Frustum3D,
    frame: Frame,
    lod: LodView,
    pyramid?: OcclusionPyramid,
  ) => boolean;
  render: (variant: GraphicsRenderVariant, globals: GlobalBindings) => void;
  /** Cached clip info, `undefined` if there are no clip objects */
  getClipInfo: () => ClipInfo | undefined;
  /** The bounding spheres of the objects drawn together (merged render objects), unscaled */
  getSegmentSpheres?: () => readonly Sphere3D[];
  getByteCount: () => number;
  getProgram: (variant: GraphicsRenderVariant) => GraphicsProgram;
  setTransparency: (transparency: Transparency) => void;
  /** Incremented when the programs change other than through the values (`setTransparency`) */
  readonly programVersion: number;
  /**
   * Advances when what decides the draw calls besides the values and the program changes: the ranges of the instances
   * and levels of detail drawn (`cull`, `uncull`, `cullSimple`), the level uniforms; equal versions mean the same calls
   * (what a recorded pass replays). Compared by those writing the ranges, a number to read
   */
  readonly drawVersion: number;
  update: () => void;
  dispose: () => void;
}

const levelIds = new WeakMap<MultiDrawBaseData, number>();
let nextLevelId = 1;

/**
 * A key written in place (`Renderable.appendDrawKey`): its first `length` values; `values` keeps its storage when
 * written again from the start (`length` set to 0), which a plain array cleared to be pushed to again does not
 */
export interface DrawKey {
  readonly values: number[];
  length: number;
}

export function DrawKey(): DrawKey {
  return { values: [], length: 0 };
}

export namespace DrawKey {
  export function add(key: DrawKey, v: number) {
    key.values[key.length++] = v;
  }

  export function areEqual(a: DrawKey, b: DrawKey) {
    if (a.length !== b.length) return false;
    const x = a.values;
    const y = b.values;
    for (let i = 0, il = a.length; i < il; ++i) {
      if (x[i] !== y[i]) return false;
    }
    return true;
  }
}

/**
 * A version of what `append` writes (`Renderable.drawVersion`): it advances when a call finds the key changed from
 * the last call's (after writing the ranges)
 */
export function createDrawVersion(append: (out: DrawKey) => void) {
  let last = DrawKey();
  let next = DrawKey();
  let version = 0;
  return () => {
    next.length = 0;
    append(next);
    if (!DrawKey.areEqual(next, last)) {
      const t = last;
      last = next;
      next = t;
      version += 1;
    }
    return version;
  };
}

const indirectBufferIds = new WeakMap<StorageBuffer, number>();
let nextIndirectBufferId = 0;
function getIndirectBufferId(buffer: StorageBuffer) {
  let id = indirectBufferIds.get(buffer);
  if (id === undefined) indirectBufferIds.set(buffer, (id = nextIndirectBufferId++));
  return id;
}

/** Appends what decides the draws of `list` (`Renderable.drawVersion`) */
export function appendMdbDrawKey(list: readonly MultiDrawBaseData[], out: DrawKey) {
  DrawKey.add(out, list.length);
  for (const l of list) {
    let id = levelIds.get(l);
    if (id === undefined) levelIds.set(l, (id = nextLevelId++));
    // a replaced level has a draw block of its own
    DrawKey.add(out, id);
    DrawKey.add(out, l.count);
    const { counts, firsts, offsets, instanceCounts, baseVertices, baseInstances, indirect } = l;
    if (indirect) {
      // the arguments are written on the GPU: the draws are the same while the buffer and offset are
      DrawKey.add(out, getIndirectBufferId(indirect.buffer));
      DrawKey.add(out, indirect.byteOffset);
    }
    for (let i = 0, il = indirect ? 0 : l.count; i < il; ++i) {
      DrawKey.add(out, counts[i]);
      DrawKey.add(out, firsts[i]);
      DrawKey.add(out, offsets[i]);
      DrawKey.add(out, instanceCounts[i]);
      DrawKey.add(out, baseVertices[i]);
      DrawKey.add(out, baseInstances[i]);
    }
    // written into the level's draw block; the cells may be refreshed with equal values, so their values count
    for (const [, cell] of l.uniforms) {
      const v = cell.ref.value as unknown;
      if (typeof v === 'number') DrawKey.add(out, v);
      else if (ArrayBuffer.isView(v) || Array.isArray(v)) {
        const a = v as ArrayLike<number>;
        for (let i = 0; i < a.length; ++i) DrawKey.add(out, a[i]);
      } else DrawKey.add(out, NaN);
    }
  }
}

function getMdbData(cellCount: number, mdbData?: MultiDrawBaseData): MultiDrawBaseData {
  if (mdbData && mdbData.instanceCounts.length >= cellCount) {
    return mdbData;
  } else {
    return {
      firsts: new Int32Array(cellCount),
      counts: new Int32Array(cellCount),
      offsets: new Int32Array(cellCount),
      instanceCounts: new Int32Array(cellCount),
      baseVertices: new Int32Array(cellCount),
      baseInstances: new Uint32Array(cellCount),
      count: 0,
      uniforms: [],
    };
  }
}

//

/**
 * Levels of detail: per level its range and fade (`overlap`), the vertices it draws and its spheres' scale. The range
 * in Å of distance from the camera plane, or with `pixelSize` above 0 (the smallest sphere drawn, in CSS pixels) in Å a
 * CSS pixel covers (the last level's maximum `Infinity`),
 * converted to distances per view (`LodView`, `screen-space-lod.md`)
 */
export type LodLevelsValue = [
  minDistance: number,
  maxDistance: number,
  overlap: number,
  count: number,
  scale: number,
][] & {
  readonly pixelSize?: number;
};

/** The largest distance of a level converted from pixel sizes (`LodView`), for an open end */
const MaxLodDistance = 1e9;

/**
 * What converts levels by size on screen (`LodLevelsValue.pixelSize`) into distances: the Å a CSS pixel covers per Å of
 * view depth (perspective), or everywhere (orthographic), from the projection and the viewport's height in CSS pixels
 * (`screen-space-lod.md`)
 */
export interface LodView {
  /** Å per CSS pixel per Å of view depth, 0 for orthographic */
  depthScale: number;
  /** Å per CSS pixel everywhere, for orthographic (0 for perspective) */
  pixelSize: number;
  /**
   * The renderer's culling view (the frustum, the camera plane and these scales), a new value where any of them changed:
   * culls of the same view and data are skipped (`createCullView`); 0 for none known (never skipped)
   */
  version: number;
}

export function LodView(): LodView {
  return { depthScale: 0, pixelSize: 0, version: 0 };
}

export namespace LodView {
  /** `2 / P[1][1]` is `2 tan(fov / 2)` (perspective) or the view's height (orthographic), `height` in CSS pixels */
  export function set(out: LodView, projection: Mat4, height: number, orthographic: boolean) {
    const s = height > 0 ? 2 / (projection[5] * height) : 0;
    out.depthScale = orthographic ? 0 : s;
    out.pixelSize = orthographic ? s : 0;
    return out;
  }

  /**
   * Level `j`'s range and fade in distances, into `out` (min, max, overlap): as they are without `pixelSize`; with it the
   * pixel sizes divided by the depth scale, or (orthographic, or no view yet) all distances for the level of the view's
   * pixel size and none for the others
   */
  export function getRange(levels: LodLevelsValue, j: number, view: LodView, out: number[]) {
    const [min, max, overlap] = levels[j];
    if (!((levels.pixelSize ?? 0) > 0)) {
      out[0] = min;
      out[1] = max;
      out[2] = overlap;
    } else if (view.depthScale > 0) {
      const k = view.depthScale;
      out[0] = min / k;
      out[1] = Math.min(max / k, MaxLodDistance);
      out[2] = overlap / k;
    } else {
      const drawn = min <= view.pixelSize && view.pixelSize < max;
      out[0] = drawn ? 0 : -2;
      out[1] = drawn ? MaxLodDistance : -1;
      out[2] = 0;
    }
    return out;
  }
}

/** The minimal set of values needed to cull a (segment of a) render item. */
export type CullValues = {
  readonly drawCount: ValueCell<number>;
  readonly instanceCount: ValueCell<number>;
  readonly uLod: ValueCell<Vec4>;
  readonly instanceGrid: ValueCell<InstanceGrid>;
  readonly boundingSphere: ValueCell<Sphere3D>;
  readonly lodLevels?: ValueCell<unknown>;
};

export type CullSegment = {
  /** first vertex of the segment, for non-indexed draws */
  first: number;
  /** byte offset into the element buffer, for indexed draws */
  offset: number;
  /** first instance of the segment in the instanced attributes */
  instanceBase: number;
};

export const ZeroCullSegment: CullSegment = { first: 0, offset: 0, instanceBase: 0 };

/**
 * A multi-draw-base-data list (one entry per lod level) that can be filled
 * from multiple segments sharing one render item.
 */
export interface SegmentedMdbList {
  readonly list: MultiDrawBaseData[];
  /** resize and reset the list and refresh the per-level uLod (and, if segmented, uLodLevel) uniforms */
  prepare: (
    lodCell: ValueCell<unknown> | undefined,
    capacity: number,
    segmented: boolean,
    lod: LodView,
  ) => LodLevelsValue | undefined;
  /** append a draw entry, extending the previous entry when contiguous */
  append: (
    l: MultiDrawBaseData,
    first: number,
    offset: number,
    count: number,
    instanceCount: number,
    baseInstance: number,
  ) => void;
  /**
   * Cull a segment against its instance grid and append draw entries, the levels' ranges those of the last `prepare`
   * (`lodLevels` gives the segment's vertex counts per level).
   * Returns false without appending if the segment has no usable grid.
   * With `wholeSegmentFallback` gridless segments are instead culled as a
   * whole using their bounding sphere.
   */
  cullSegment: (
    values: CullValues,
    segment: CullSegment,
    lodLevels: LodLevelsValue | undefined,
    cameraPlane: Plane3D,
    frustum: Frustum3D,
    isOccluded: ((s: Sphere3D) => boolean) | null,
    stats: GraphicsStats,
    wholeSegmentFallback: boolean,
  ) => boolean;
  /** append draw entries for the whole segment without culling */
  fullSegment: (values: CullValues, segment: CullSegment, lodLevels: LodLevelsValue | undefined) => void;
  /** append the draw entry for the lod level matching the given distance (the ranges of the last `prepare`) */
  cullSimpleSegment: (
    values: CullValues,
    segment: CullSegment,
    lodLevels: LodLevelsValue,
    d: number,
    radius: number,
    scale: number,
  ) => void;
}

/**
 * The frustum's planes and the camera plane as coefficients (`[nx, ny, nz, constant]` each, the camera plane last), set
 * per cull: the sphere tests run per batch and cell of the instance grids, reading the spheres straight from their arrays
 * (through `Frustum3D`, `Plane3D` and `Sphere3D` they cost about 90% of a rotating frame of 100k instances)
 */
const cullPlanes = new Float64Array(28);

function setCullPlanes(frustum: Frustum3D, cameraPlane: Plane3D) {
  for (let i = 0; i < 6; ++i) {
    const { normal, constant } = frustum[i as 0];
    cullPlanes[i * 4] = normal[0];
    cullPlanes[i * 4 + 1] = normal[1];
    cullPlanes[i * 4 + 2] = normal[2];
    cullPlanes[i * 4 + 3] = constant;
  }
  cullPlanes[24] = cameraPlane.normal[0];
  cullPlanes[25] = cameraPlane.normal[1];
  cullPlanes[26] = cameraPlane.normal[2];
  cullPlanes[27] = cameraPlane.constant;
}

/** As `Frustum3D.intersectsSphere3D` with the planes of `setCullPlanes` */
function intersectsCullFrustum(x: number, y: number, z: number, r: number) {
  const p = cullPlanes;
  for (let i = 0; i < 24; i += 4) {
    if (p[i] * x + p[i + 1] * y + p[i + 2] * z + p[i + 3] < -r) return false;
  }
  return true;
}

/** As `Frustum3D.containsSphere3D` with the planes of `setCullPlanes` */
function containedInCullFrustum(x: number, y: number, z: number, r: number) {
  const p = cullPlanes;
  for (let i = 0; i < 24; i += 4) {
    if (p[i] * x + p[i + 1] * y + p[i + 2] * z + p[i + 3] < r) return false;
  }
  return true;
}

/** As `Plane3D.distanceToPoint` of the camera plane of `setCullPlanes` */
function cullCameraDistance(x: number, y: number, z: number) {
  const p = cullPlanes;
  return p[24] * x + p[25] * y + p[26] * z + p[27];
}

export function createSegmentedMdbList(): SegmentedMdbList {
  const list: MultiDrawBaseData[] = [];
  let lodLevelsVersion = -1;
  // the view the level uniforms were converted with (levels by pixel size)
  let lodDepthScale = -1;
  let lodPixelSize = -1;
  const range: number[] = [0, 0, 0];

  const s = Sphere3D();

  // flattened lod tuples, avoids nested array loads in the hot cell loop; the ranges in distances (`prepare`)
  const lodMin: number[] = [];
  const lodMax: number[] = [];
  const lodOverlap: number[] = [];
  const lodCount: number[] = [];

  function prepare(
    lodCell: ValueCell<unknown> | undefined,
    capacity: number,
    segmented: boolean,
    lod: LodView,
  ): LodLevelsValue | undefined {
    const lodLevels = lodCell?.ref.value as LodLevelsValue | undefined;
    const hasLodLevels = !!lodLevels && lodLevels.length > 0;
    const levelCount = hasLodLevels ? lodLevels.length : 1;
    // segmented (merged) items also carry a per-instance `aSegmentLod`
    // (up to 4 levels) - members' scales can differ, unlike min/max/overlap
    const uniformCount = segmented ? 2 : 1;

    list.length = levelCount;
    const viewChanged =
      (lodLevels?.pixelSize ?? 0) > 0 && (lod.depthScale !== lodDepthScale || lod.pixelSize !== lodPixelSize);
    const uniformsNeedUpdate = hasLodLevels && (lodCell!.ref.version !== lodLevelsVersion || viewChanged);
    for (let j = 0; j < levelCount; ++j) {
      list[j] = getMdbData(capacity, list[j]);
      list[j].count = 0;
      if (hasLodLevels) {
        LodView.getRange(lodLevels, j, lod, range);
        lodMin[j] = range[0];
        lodMax[j] = range[1];
        lodOverlap[j] = range[2];
        const l = list[j];
        let created = false;
        if (l.uniforms.length !== uniformCount) {
          l.uniforms.length = uniformCount;
          l.uniforms[0] = ['uLod', ValueCell.create(Vec4())];
          if (segmented) l.uniforms[1] = ['uLodLevel', ValueCell.create(Math.min(j, 3))];
          created = true;
        }
        if (uniformsNeedUpdate || created) {
          ValueCell.update(
            l.uniforms[0][1],
            Vec4.set(l.uniforms[0][1].ref.value as Vec4, lodMin[j], lodMax[j], lodOverlap[j], lodLevels[j][4]),
          );
        }
      } else {
        list[j].uniforms.length = 0;
      }
    }
    if (uniformsNeedUpdate) {
      lodLevelsVersion = lodCell!.ref.version;
      lodDepthScale = lod.depthScale;
      lodPixelSize = lod.pixelSize;
    }

    return hasLodLevels ? lodLevels : undefined;
  }

  function emit(
    l: MultiDrawBaseData,
    first: number,
    offset: number,
    count: number,
    instanceCount: number,
    baseInstance: number,
  ) {
    const o = l.count;
    if (
      o > 0 &&
      l.firsts[o - 1] === first &&
      l.offsets[o - 1] === offset &&
      l.counts[o - 1] === count &&
      l.baseInstances[o - 1] + l.instanceCounts[o - 1] === baseInstance
    ) {
      l.instanceCounts[o - 1] += instanceCount;
    } else {
      l.firsts[o] = first;
      l.offsets[o] = offset;
      l.counts[o] = count;
      l.instanceCounts[o] = instanceCount;
      l.baseInstances[o] = baseInstance;
      l.count += 1;
    }
  }

  function cullSegment(
    values: CullValues,
    segment: CullSegment,
    lodLevels: LodLevelsValue | undefined,
    cameraPlane: Plane3D,
    frustum: Frustum3D,
    isOccluded: ((s: Sphere3D) => boolean) | null,
    stats: GraphicsStats,
    wholeSegmentFallback: boolean,
  ): boolean {
    const drawCount = values.drawCount.ref.value;
    const instanceCount = values.instanceCount.ref.value;
    if (drawCount === 0 || instanceCount === 0) return false;

    const { first, offset, instanceBase } = segment;
    const [minDistance, maxDistance] = values.uLod.ref.value;
    const hasLod = minDistance !== 0 || maxDistance !== 0;
    const grid = values.instanceGrid.ref.value;

    const levelCount = lodLevels ? lodLevels.length : 0;
    for (let j = 0; j < levelCount; ++j) lodCount[j] = lodLevels![j][3];

    if (grid.cellSize <= 1) {
      if (!wholeSegmentFallback) return false;

      // no instance grid for this segment, cull it as a whole
      const b = values.boundingSphere.ref.value;
      const d = p3distanceToPoint(cameraPlane, b.center);
      if (hasLod && (d + b.radius < minDistance || d - b.radius > maxDistance)) {
        if (isTimingMode) stats.culled.lod += instanceCount;
        return true;
      }
      if (!f3intersectsSphere3D(frustum, b)) {
        if (isTimingMode) stats.culled.frustum += instanceCount;
        return true;
      }
      if (isOccluded !== null && isOccluded(b)) {
        if (isTimingMode) stats.culled.occlusion += instanceCount;
        return true;
      }
      if (lodLevels) {
        // no per-cell distance to decimate by here, so draw full detail for
        // whichever level(s) match (keeps the level's uLod min/max/overlap
        // fade correct, matches the plain render path's full-detail draw for
        // ungridded objects, see the `cellSize <= 1` bail in createRenderable's cull)
        for (let j = 0, jl = lodLevels.length; j < jl; ++j) {
          if (d + b.radius < lodMin[j] || d - b.radius > lodMax[j]) continue;
          emit(list[j], first, offset, drawCount, instanceCount, instanceBase);
        }
      } else {
        emit(list[0], first, offset, drawCount, instanceCount, instanceBase);
      }
      return true;
    }

    const { cellOffsets, cellSpheres, batchOffsets, batchSpheres, batchCount, batchSize } = grid;
    const checkCellOccludedDistance = 2 * batchSize;
    setCullPlanes(frustum, cameraPlane);

    for (let k = 0; k < batchCount; ++k) {
      // batchCell is serial after the grid reorder, cells of a batch are consecutive
      const cBegin = batchOffsets[k];
      const cEnd = batchOffsets[k + 1];
      const cCount = cEnd - cBegin;
      if (cCount === 0) continue;

      const bx = batchSpheres[k * 4];
      const by = batchSpheres[k * 4 + 1];
      const bz = batchSpheres[k * 4 + 2];
      const br = batchSpheres[k * 4 + 3];
      const bd = cullCameraDistance(bx, by, bz);
      if (hasLod && (bd + br < minDistance || bd - br > maxDistance)) {
        if (isTimingMode) {
          stats.culled.lod += cellOffsets[cEnd] - cellOffsets[cBegin];
        }
        continue;
      }
      if (!intersectsCullFrustum(bx, by, bz, br)) {
        if (isTimingMode) {
          stats.culled.frustum += cellOffsets[cEnd] - cellOffsets[cBegin];
        }
        continue;
      }
      if (isOccluded !== null && isOccluded(s3fromArray(s, batchSpheres, k * 4))) {
        if (isTimingMode) {
          stats.culled.occlusion += cellOffsets[cEnd] - cellOffsets[cBegin];
        }
        continue;
      }

      if (cCount === 1) {
        // single-cell batch: cell sphere equals batch sphere by construction
        // (see calcTopGrid), re-testing it per-cell would be redundant
        const begin = cellOffsets[cBegin];
        const count = cellOffsets[cBegin + 1] - begin;
        if (count === 0) continue;

        if (lodLevels) {
          const dMax = bd + br;
          const dMin = bd - br;
          for (let j = 0; j < levelCount; ++j) {
            if (dMax < lodMin[j] || dMin > lodMax[j]) continue;
            emit(list[j], first, offset, lodCount[j], count, instanceBase + begin);
          }
        } else {
          emit(list[0], first, offset, drawCount, count, instanceBase + begin);
        }
        continue;
      }

      // cells of a fully contained batch don't need frustum tests
      const batchContained = containedInCullFrustum(bx, by, bz, br);

      // a contained batch drawn by its cells as a whole: its instances at once, without testing each cell (the
      // cells' instances are consecutive, `emit` joins them all the same); where no cell tests occlusion and the
      // batch lies within the distance range and one level alone (each cell lies within the batch's sphere)
      if (batchContained && (isOccluded === null || bd - br >= checkCellOccludedDistance)) {
        const dMin = bd - br;
        const dMax = bd + br;
        if (!hasLod || (dMin >= minDistance && dMax <= maxDistance)) {
          let level = lodLevels ? -1 : 0;
          if (lodLevels) {
            for (let j = 0; j < levelCount; ++j) {
              if (dMax < lodMin[j] || dMin > lodMax[j]) continue;
              // overlaps a second level, or not within the one
              if (level !== -1 || dMin < lodMin[j] || dMax > lodMax[j]) {
                level = -2;
                break;
              }
              level = j;
            }
          }
          if (level >= 0) {
            const begin = cellOffsets[cBegin];
            const count = cellOffsets[cEnd] - begin;
            if (count > 0) {
              const vertices = lodLevels ? lodCount[level] : drawCount;
              emit(list[level], first, offset, vertices, count, instanceBase + begin);
            }
            continue;
          }
        }
      }

      for (let c = cBegin; c < cEnd; ++c) {
        const begin = cellOffsets[c];
        const end = cellOffsets[c + 1];
        const count = end - begin;
        if (count === 0) continue;

        const x = cellSpheres[c * 4];
        const y = cellSpheres[c * 4 + 1];
        const z = cellSpheres[c * 4 + 2];
        const r = cellSpheres[c * 4 + 3];
        const d = cullCameraDistance(x, y, z);
        if (hasLod && (d + r < minDistance || d - r > maxDistance)) {
          if (isTimingMode) stats.culled.lod += count;
          continue;
        }
        if (!batchContained && !intersectsCullFrustum(x, y, z, r)) {
          if (isTimingMode) stats.culled.frustum += count;
          continue;
        }
        if (
          isOccluded !== null &&
          d - r < checkCellOccludedDistance &&
          isOccluded(s3fromArray(s, cellSpheres, c * 4))
        ) {
          if (isTimingMode) stats.culled.occlusion += count;
          continue;
        }

        if (lodLevels) {
          const dMax = d + r;
          const dMin = d - r;
          for (let j = 0; j < levelCount; ++j) {
            if (dMax < lodMin[j] || dMin > lodMax[j]) continue;
            emit(list[j], first, offset, lodCount[j], count, instanceBase + begin);
          }
        } else {
          emit(list[0], first, offset, drawCount, count, instanceBase + begin);
        }
      }
    }
    return true;
  }

  function fullSegment(values: CullValues, segment: CullSegment, lodLevels: LodLevelsValue | undefined) {
    const drawCount = values.drawCount.ref.value;
    const instanceCount = values.instanceCount.ref.value;
    if (drawCount === 0 || instanceCount === 0) return;

    if (lodLevels) {
      for (let j = 0, jl = lodLevels.length; j < jl; ++j) {
        emit(list[j], segment.first, segment.offset, lodLevels[j][3], instanceCount, segment.instanceBase);
      }
    } else {
      emit(list[0], segment.first, segment.offset, drawCount, instanceCount, segment.instanceBase);
    }
  }

  function cullSimpleSegment(
    values: CullValues,
    segment: CullSegment,
    lodLevels: LodLevelsValue,
    d: number,
    radius: number,
    scale: number,
  ) {
    if (values.drawCount.ref.value === 0 || values.instanceCount.ref.value === 0) return;

    for (let j = 0, jl = lodLevels.length; j < jl; ++j) {
      if (d + radius < lodMax[j] * scale) {
        emit(
          list[j],
          segment.first,
          segment.offset,
          lodLevels[j][3],
          values.instanceCount.ref.value,
          segment.instanceBase,
        );
        break;
      }
    }
  }

  return { list, prepare, append: emit, cullSegment, fullSegment, cullSimpleSegment };
}

type GraphicsRenderableValues = RenderableValues & BaseValues;

/**
 * What a cull read, to skip the next with the same (still frames redo every object's culling otherwise): the view's
 * version (`LodView.version`, the renderer compares the view once per update) and the versions of the culled values'
 * cells. Not for an occlusion test (its levels change between frames) nor in timing mode (the culled counts are per
 * frame); `invalidate` where the lists changed otherwise (`uncull`, `cullSimple`)
 */
export function createCullView() {
  let viewVersion = 0;
  const versions: number[] = [];
  let valid = false;
  return {
    /** whether a cull with these would cull as the last did */
    isSame(lod: LodView, isOccluded: unknown, cells: readonly (ValueCell<unknown> | undefined)[]) {
      if (!valid || lod.version === 0 || lod.version !== viewVersion || isOccluded !== null || isTimingMode)
        return false;
      for (let i = 0, il = cells.length; i < il; ++i) if ((cells[i]?.ref.version ?? -1) !== versions[i]) return false;
      return true;
    },
    set(lod: LodView, cells: readonly (ValueCell<unknown> | undefined)[]) {
      viewVersion = lod.version;
      versions.length = cells.length;
      for (let i = 0, il = cells.length; i < il; ++i) versions[i] = cells[i]?.ref.version ?? -1;
      valid = true;
    },
    invalidate() {
      valid = false;
    },
  };
}

/** The cells of `values` culling reads (`createCullView`) */
export function getCullCells(values: CullValues): (ValueCell<unknown> | undefined)[] {
  return [
    values.drawCount,
    values.instanceCount,
    values.uLod,
    values.instanceGrid,
    values.boundingSphere,
    values.lodLevels,
  ];
}

export function createRenderable<T extends GraphicsRenderableValues>(
  renderItem: GraphicsRenderItem,
  values: T,
  state: RenderableState,
): Renderable<T> {
  const id = getNextRenderableId();

  const mdb = createSegmentedMdbList();
  const clipInfo = createClipInfoCache(values);
  let cullEnabled = false;
  let lastCullFrame: Frame | undefined;
  let programVersion = 0;
  const getDrawVersion = createDrawVersion((out) => {
    DrawKey.add(out, cullEnabled ? 1 : 0);
    if (cullEnabled) appendMdbDrawKey(mdb.list, out);
  });
  let drawVersion = getDrawVersion();
  // the cells culling reads every frame, resolved once: values objects have too many keys for fast named reads (the
  // cells stay, updates replace their values)
  const cullValues: CullValues = {
    drawCount: values.drawCount,
    instanceCount: values.instanceCount,
    uLod: values.uLod,
    instanceGrid: values.instanceGrid,
    boundingSphere: values.boundingSphere,
    lodLevels: values.lodLevels,
  };

  const cullView = createCullView();
  const cullCells = getCullCells(cullValues);

  const cull = (
    cameraPlane: Plane3D,
    frustum: Frustum3D,
    isOccluded: ((s: Sphere3D) => boolean) | null,
    stats: GraphicsStats,
    frame: Frame,
    lod: LodView,
  ) => {
    cullEnabled = false;
    cullView.set(lod, cullCells);

    if (cullValues.drawCount.ref.value === 0) return;
    if (cullValues.instanceCount.ref.value === 0) return;
    const grid = cullValues.instanceGrid.ref.value;
    if (grid.cellSize <= 1) return;

    const lodLevels = mdb.prepare(cullValues.lodLevels, grid.cellCount, false, lod);
    mdb.cullSegment(cullValues, ZeroCullSegment, lodLevels, cameraPlane, frustum, isOccluded, stats, false);

    cullEnabled = true;
    lastCullFrame = frame;
  };

  const cullSimple = (d: number, radius: number, scale: number, lod: LodView) => {
    cullEnabled = false;

    const lodLevels = mdb.prepare(
      cullValues.lodLevels,
      Math.max(1, cullValues.instanceGrid.ref.value.cellCount),
      false,
      lod,
    );
    if (!lodLevels) return;

    mdb.cullSimpleSegment(cullValues, ZeroCullSegment, lodLevels, d, radius, scale);

    cullEnabled = true;
  };

  return {
    id,
    materialId: renderItem.materialId,
    renderItem,
    values,
    state,

    cull: (
      cameraPlane: Plane3D,
      frustum: Frustum3D,
      isOccluded: ((s: Sphere3D) => boolean) | null,
      stats: GraphicsStats,
      frame: Frame,
      lod: LodView,
    ) => {
      // skip recomputation if this renderable was already culled for the current frame, or for the same view and data
      if (frame === lastCullFrame) {
        stats.cacheHits.cull++;
        return;
      }
      if (cullView.isSame(lod, isOccluded, cullCells)) {
        lastCullFrame = frame;
        stats.cacheHits.cull++;
        return;
      }
      cull(cameraPlane, frustum, isOccluded, stats, frame, lod);
      drawVersion = getDrawVersion();
    },
    uncull: () => {
      lastCullFrame = undefined;
      cullView.invalidate();
      if (!cullEnabled) return;
      cullEnabled = false;
      drawVersion = getDrawVersion();
    },
    cullSimple: (d: number, radius: number, scale: number, lod: LodView) => {
      cullView.invalidate();
      cullSimple(d, radius, scale, lod);
      drawVersion = getDrawVersion();
    },
    render: (variant: GraphicsRenderVariant, globals: GlobalBindings) => {
      renderItem.render(variant, globals, cullEnabled ? mdb.list : undefined);
    },
    getClipInfo: clipInfo.get,
    getByteCount: () => renderItem.getByteCount(),
    getProgram: (variant: GraphicsRenderVariant) => renderItem.getProgram(variant),
    setTransparency: (transparency: Transparency) => {
      renderItem.setTransparency(transparency);
      programVersion += 1;
    },
    get programVersion() {
      return programVersion;
    },
    get drawVersion() {
      return drawVersion;
    },
    update: () => {
      renderItem.update();
      lastCullFrame = undefined;
    },
    dispose: () => renderItem.destroy(),
  };
}

export type GraphicsRenderable = Renderable<GraphicsRenderableValues>;

//

export type ClipInfo = {
  objects: Clip.Objects;
  /** Planes of the plane-type clip objects, see `Clip.getPlane`, only for the pixel variant */
  readonly planes: Plane3D[];
  /** Indices of clip objects whose surface can intersect the bounding sphere, only for the pixel variant */
  readonly capIndices: number[];
};

export interface ClipInfoCache {
  /** Clip info, recomputed when any of its inputs changed, `undefined` if there are no clip objects */
  readonly get: () => ClipInfo | undefined;
}

function updateClipInfo(info: ClipInfo, values: BaseValues) {
  const { planes, capIndices } = info;
  const objects = values.clipObjects.ref.value;
  info.objects = objects;

  capIndices.length = 0;
  if (values.dClipVariant.ref.value !== 'pixel') return;

  const boundingSphere = values.boundingSphere.ref.value;
  for (let i = 0; i < objects.count; ++i) {
    if (objects.type[i] === Clip.Type.plane) Clip.getPlane((planes[i] ??= Plane3D()), objects, i);
    if (Clip.canIntersectSphere(objects, i, boundingSphere)) capIndices.push(i);
  }
}

export function createClipInfoCache(values: BaseValues): ClipInfoCache {
  const cells: ValueCell<unknown>[] = [
    values.uClipObjectCount,
    values.dClipVariant,
    values.boundingSphere,
    values.clipObjects,
  ];
  const versions = cells.map(() => -1);

  const info: ClipInfo = {
    objects: values.clipObjects.ref.value,
    planes: [],
    capIndices: [],
  };

  function get(): ClipInfo | undefined {
    if (values.uClipObjectCount.ref.value === 0) return;

    let changed = false;
    for (let i = 0, il = cells.length; i < il; ++i) {
      const version = cells[i].ref.version;
      if (version !== versions[i]) {
        versions[i] = version;
        changed = true;
      }
    }
    if (changed) updateClipInfo(info, values);
    return info;
  }

  return { get };
}

//

export interface ComputeRenderable<T extends RenderableValues> {
  readonly id: number;
  readonly values: T;
  /** For pipeline requests (`resources.preparePipelines`) */
  readonly renderItem: ComputeRenderItem;

  render: () => void;
  update: () => void;
  dispose: () => void;
}

export function createComputeRenderable<T extends Values<RenderableSchema>>(
  renderItem: ComputeRenderItem,
  values: T,
): ComputeRenderable<T> {
  return {
    id: getNextRenderableId(),
    values,
    renderItem,

    render: () => {
      renderItem.getProgram('compute').finalize(true);
      renderItem.render('compute');
    },
    update: () => renderItem.update(),
    dispose: () => renderItem.destroy(),
  };
}
