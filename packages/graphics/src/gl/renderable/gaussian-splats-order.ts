/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '@molstar/core/math/linear-algebra';
import { countingSortIndices } from '@molstar/core/data/util/sort';
import { ValueCell } from '@molstar/core/util';
import { type DataArray, createDataArray } from './util.js';
import type { ComputeContext, StorageBuffer } from '../compute.js';
import type { SplatSorter, SplatSortJob } from '../compute/splat-sort.js';
import type { GaussianSplatsValues } from './gaussian-splats.js';

/**
 * The draw order of gaussian splats (`tSplatOrder`): per instance a block of the splats' indices, back to front for the
 * view it was sorted for; a `u32` an index as the four bytes of an `rgba8` element (little-endian), read by the vertex
 * shader for its draw slot (`aInstance * uVertexCount + slot`)
 */
export function createSplatOrder(count: number, instanceCount: number, array?: Uint8Array): DataArray<Uint8Array> {
  const data = createDataArray(count * instanceCount, 4, Uint8Array, array);
  const words = getSplatOrderWords(data);
  for (let k = 0; k < instanceCount; ++k) {
    for (let i = 0, o = k * count; i < count; ++i) words[o + i] = i;
  }
  return data;
}

/** The indices of `data` (`createSplatOrder`) as `u32` */
export function getSplatOrderWords(data: DataArray<Uint8Array>) {
  const { array, count } = data;
  return new Uint32Array(array.buffer, array.byteOffset, count);
}

/** Steps of the quantized depth: a 16-bit key */
const KeySteps = 65535;

/**
 * Writes `count` splat indices into `out` from `offset`, back to front along the view row `(a, b, c)` (a splat's view
 * z, up to a constant, is `a x + b y + c z` of its center: the farthest the most negative, first); centers from
 * `positionGroup` (`PositionGroupStride` floats a splat); `keys` scratch of at least `count`
 */
export function sortSplats(
  out: Uint32Array,
  offset: number,
  positionGroup: Float32Array,
  count: number,
  a: number,
  b: number,
  c: number,
  keys: Uint32Array,
) {
  if (count === 0) return;
  let min = Infinity;
  let max = -Infinity;
  // the depths first, into `keys` as floats' bits would lose the order: kept in a view of its buffer; their bounds of
  // the depths as stored (float32: a bound of the float64 ones would key a rounded depth past the keys' range)
  const depths = new Float32Array(keys.buffer, keys.byteOffset, count);
  for (let i = 0, o = 0; i < count; ++i, o += 4) {
    depths[i] = a * positionGroup[o] + b * positionGroup[o + 1] + c * positionGroup[o + 2];
    const z = depths[i];
    if (z < min) min = z;
    if (z > max) max = z;
  }
  const scale = max > min ? KeySteps / (max - min) : 0;
  // in place: a key written after its depth is read (the same index); within the keys' range (a NaN center's 0)
  for (let i = 0; i < count; ++i) {
    const key = ((depths[i] - min) * scale) | 0;
    keys[i] = key < 0 ? 0 : key > KeySteps ? KeySteps : key;
  }
  // one pass reading the keys in order (1 M splats 8-10 ms against 22.5 of `radixSortIndices` in Node, 5 M 42-47 against 188)
  countingSortIndices(keys, count, KeySteps, out, offset);
}

/** Below this cosine between the view rows' directions an instance is sorted again (about 1.8°) */
const ResortCosine = 0.9995;

interface SortState {
  /** per instance the unit direction of the view row last sorted for (NaN: never) */
  directions: Float32Array;
  keys: Uint32Array;
  /** the versions sorted for: centers, transforms */
  positionVersion: number;
  transformVersion: number;
  /** sorted on the GPU last (into `sSplatOrder`), else on the CPU (into `tSplatOrder`) */
  gpu: boolean;
  /**
   * of the GPU sort: the centers (`vec4` a splat) and what they were made from (the CPU's centers' version, or SOG's
   * means: decoded on the GPU), the order, those replaced (destroyed a sort later)
   */
  positions?: StorageBuffer;
  positionsVersion: number | object;
  order?: StorageBuffer;
  retired: StorageBuffer[];
}

const SortStates = new WeakMap<object, SortState>();

const tmpInstance = Mat4();
const tmpObject = Mat4();
const tmpView = Mat4();

/**
 * Forgets the GPU's order of `values` where it is of another context than `graphics` (a backend switch, a restore):
 * the order unset (`dSplatGpuOrder` off, `tSplatOrder` read till sorted again), the sort's state and buffers (the old
 * device's, not destroyed) dropped, so the next sort is a whole one on `graphics`; true where it did
 */
export function releaseForeignSplatOrder(values: GaussianSplatsValues, compute: ComputeContext | undefined): boolean {
  const order = values.sSplatOrder.ref.value as StorageBuffer | undefined;
  if (!order || compute?.ownsBuffer(order)) return false;
  ValueCell.update(values.dSplatGpuOrder, false);
  ValueCell.update(values.sSplatOrder, undefined as unknown as StorageBuffer);
  SortStates.delete(values.tSplatOrder);
  return true;
}

/**
 * Destroys the GPU sort's buffers of `values` (its renderable disposed: drawn by no command still to be submitted) and
 * forgets its sort
 */
export function disposeSplatOrder(values: GaussianSplatsValues) {
  const state = SortStates.get(values.tSplatOrder);
  if (!state) return;
  for (const b of state.retired) b.destroy();
  state.positions?.destroy();
  state.order?.destroy();
  SortStates.delete(values.tSplatOrder);
}

/** The GPU's sort of splats (WebGPU), see `updateSplatOrder` */
export interface SplatOrderGpu {
  readonly compute: ComputeContext;
  readonly sorter: SplatSorter;
}

/**
 * Sorts the splats of `values` back to front for the model-view `modelView` (the camera's view and the scene's model;
 * a uniform scale changes no order), each instance whose view row turned since its last sort (z sorting: moving
 * without turning keeps the order): on the GPU with `gpu` (into `sSplatOrder`, `dSplatGpuOrder` set), else on the CPU
 * (into `tSplatOrder`); true where the values changed (the renderable to be updated)
 */
export function updateSplatOrder(values: GaussianSplatsValues, modelView: Mat4, gpu?: SplatOrderGpu): boolean {
  const pg = values.splatCenters.ref.value;
  const count = values.uVertexCount.ref.value;
  const instanceCount = values.instanceCount.ref.value;
  if (count < 2) return false;
  const orderCell = values.tSplatOrder;
  let data = orderCell.ref.value;
  let resized = false;
  if (data.count !== count * instanceCount) {
    data = createSplatOrder(count, instanceCount, data.array);
    resized = true;
  }

  let state = SortStates.get(orderCell);
  if (!state || state.directions.length !== instanceCount * 3 || resized) {
    state = {
      directions: new Float32Array(instanceCount * 3).fill(Number.NaN),
      keys: state && state.keys.length >= count ? state.keys : new Uint32Array(count),
      positionVersion: -1,
      transformVersion: -1,
      gpu: !!gpu,
      positions: state?.positions,
      positionsVersion: state?.positionsVersion ?? -1,
      order: state?.order,
      retired: state?.retired ?? [],
    };
    SortStates.set(orderCell, state);
  }
  const all =
    state.positionVersion !== values.splatCenters.ref.version ||
    state.transformVersion !== values.aTransform.ref.version ||
    state.gpu !== !!gpu;
  state.gpu = !!gpu;

  const words = getSplatOrderWords(data);
  const transforms = values.aTransform.ref.value;
  // row `r` of the transforms instance `aInstance[r]` (an instance grid stores them in its order)
  const instances = values.aInstance.ref.value;
  Mat4.mul(tmpObject, modelView, values.uObjectTransform.ref.value);
  const sphere = values.invariantBoundingSphere.ref.value;
  const jobs: SplatSortJob[] = [];
  let changed = resized;
  for (let r = 0; r < instanceCount; ++r) {
    const k = instances[r];
    Mat4.fromArray(tmpInstance, transforms, r * 16);
    const m = Mat4.mul(tmpView, tmpObject, tmpInstance);
    // the view z row: column-major, row 2 of the upper 3x3
    const a = m[2];
    const b = m[6];
    const c = m[10];
    const length = Math.sqrt(a * a + b * b + c * c) || 1;
    const d = state.directions;
    const o = k * 3;
    if (!all && d[o] * a + d[o + 1] * b + d[o + 2] * c >= ResortCosine * length) continue;
    if (gpu) {
      // the depths within the bounds: their center's, give or take the radius (the row's length its scale)
      const { center, radius } = sphere;
      const z = a * center[0] + b * center[1] + c * center[2];
      jobs.push({ offset: k * count, count, a, b, c, depthMin: z - radius * length, depthMax: z + radius * length });
    } else {
      sortSplats(words, k * count, pg.array, count, a, b, c, state.keys);
      changed = true;
    }
    d[o] = a / length;
    d[o + 1] = b / length;
    d[o + 2] = c / length;
  }
  state.positionVersion = values.splatCenters.ref.version;
  state.transformVersion = values.aTransform.ref.version;
  if (gpu) {
    // the order of the new counts kept (else the next sort finds them resized again, all sorted again: every frame)
    if (resized) ValueCell.update(orderCell, data);
    return sortOnGpu(values, state, gpu, count * instanceCount, jobs) || resized;
  }
  if (values.dSplatGpuOrder.ref.value) {
    ValueCell.update(values.dSplatGpuOrder, false);
    changed = true;
  }
  if (changed) ValueCell.update(orderCell, data);
  return changed;
}

/** The GPU's sort of `jobs` into the state's order (`length` indices), its buffers made or grown first */
function sortOnGpu(
  values: GaussianSplatsValues,
  state: SortState,
  gpu: SplatOrderGpu,
  length: number,
  jobs: SplatSortJob[],
): boolean {
  const { compute, sorter } = gpu;
  // replaced at the last sort: the commands reading them submitted since
  for (const b of state.retired) b.destroy();
  state.retired.length = 0;
  let changed = false;
  const pg = values.splatCenters.ref.value;
  const count = values.uVertexCount.ref.value;
  if (!state.positions || state.positions.byteLength < count * 16) {
    if (state.positions) state.retired.push(state.positions);
    state.positions = compute.createBuffer(count * 16, 'splat centers');
    state.positionsVersion = -1;
  }
  if (values.dSplatEncoding.ref.value === 'sog') {
    // SOG read in place: the centers decoded from its means on the GPU, not uploaded (the CPU's may come later)
    const { tSogMeansLow, tSogMeansHigh, uSogMeansMin, uSogMeansMax } = values;
    const source = tSogMeansLow.ref.value;
    if (state.positionsVersion !== source) {
      const min = uSogMeansMin.ref.value;
      const max = uSogMeansMax.ref.value;
      gpu.sorter.decodeSogCenters(source, tSogMeansHigh.ref.value, min, max, count, state.positions);
      state.positionsVersion = source;
    }
  } else if (state.positionsVersion !== values.splatCenters.ref.version) {
    state.positions.write(pg.array.subarray(0, count * 4));
    state.positionsVersion = values.splatCenters.ref.version;
  }
  if (!state.order || state.order.byteLength < length * 4) {
    // the old one may still be read by commands not yet submitted: destroyed at the next sort
    if (state.order) state.retired.push(state.order);
    state.order = compute.createBuffer(length * 4, 'splat order');
    ValueCell.update(values.sSplatOrder, state.order);
    changed = true;
  }
  if (jobs.length) sorter.sort(state.positions, state.order, jobs);
  if (!values.dSplatGpuOrder.ref.value) {
    ValueCell.update(values.dSplatGpuOrder, true);
    changed = true;
  }
  return changed;
}
