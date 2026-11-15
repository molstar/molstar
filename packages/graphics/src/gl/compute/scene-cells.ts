/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import type { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import {
  type ComputeContext,
  type ComputeKernel,
  KernelCode,
  type KernelRequest,
  type StorageBuffer,
  StorageBufferSpec,
} from '../compute.js';
import { UniformSpec } from '../renderable/schema.js';

/**
 * The scene cells of a merged object culled on the GPU (`merged-culler.ts`), built on the GPU from the instances' data
 * (`CulledData`): one grid over all instances by position, a grid per size class (`getSceneGrid`, on the CPU from the
 * members' bounds). `scene-cell-bin` (a thread per instance) writes each instance's cell and segment and grows its
 * cell's box of centers and largest radius by atomics (floats as ordered `u32`s); `scene-cell-spheres` (a thread per
 * cell) writes each cell's sphere from its box, an empty cell's radius negative. The cells are the grids' dense cells,
 * not compacted: their count known on the CPU (no read back), the empty ones a thread of the cell test each
 */
export const SceneCellBinKernel = KernelCode('scene-cell-bin', [64]);
export const SceneCellSpheresKernel = KernelCode('scene-cell-spheres', [64]);

/** u32 words of a cell's box (`cellBoxes`): its centers' min (inverted) and max as ordered bits, its largest radius */
export const SceneCellBoxWords = 7;

export const SceneCellBinSchema = {
  uInstanceCount: UniformSpec('i'),
  /** the grids' lower corner (xyz) */
  uSceneMin: UniformSpec('v4'),
  /** per class (column) its cell side's inverse and its dimensions */
  uSceneGrids: UniformSpec('m4'),
  /** per class its first cell */
  uSceneFirsts: UniformSpec('iv4'),
  /** per segment three vec4s, the first its invariant sphere (the culling's `segmentData`) */
  segmentData: StorageBufferSpec('vec4<f32>', 'read'),
  /** per segment its class (`SceneGrid.classes`) */
  cellClasses: StorageBufferSpec('u32', 'read'),
  /** the instances' data (`CulledData`) */
  data: StorageBufferSpec('vec4<f32>', 'read'),
  /** per instance its cell and segment, written */
  instanceCells: StorageBufferSpec('vec2<u32>', 'read-write'),
  /** per cell its box (`SceneCellBoxWords`), zeroed before */
  cellBoxes: StorageBufferSpec('atomic<u32>', 'read-write'),
};

export const SceneCellSpheresSchema = {
  uCellCount: UniformSpec('i'),
  cellBoxes: StorageBufferSpec('u32', 'read'),
  /** per cell its sphere: its box's, expanded by the largest radius; empty: a negative radius */
  cellSpheres: StorageBufferSpec('vec4<f32>', 'read-write'),
};

export const SceneCellKernels: readonly KernelRequest[] = [
  { code: SceneCellBinKernel, schema: SceneCellBinSchema },
  { code: SceneCellSpheresKernel, schema: SceneCellSpheresSchema },
];

/** The smallest cell side and scene extent (Å), for instances without a radius (the side follows their density) */
const MinSceneCell = 1;
/** Cell sizes a scene's grids have: a member's instances in the grid whose cells are at least their diameter */
export const SceneCellClasses = 4;
/** Cells a class's grid has at most (its side grown where the scene is larger) */
const MaxSceneCells = 1 << 22;

/** The grids of a merged object's scene cells, a grid per size class over the scene's box (`getSceneGrid`) */
export interface SceneGrid {
  /** the box's lower corner */
  readonly min: Vec4;
  /** per class (column) its cell side's inverse and its dimensions (`uSceneGrids`) */
  readonly grids: Mat4;
  /** per class its first cell (`uSceneFirsts`) */
  readonly firsts: Vec4;
  /** the cells of all classes */
  readonly cellCount: number;
  /** per segment its class */
  readonly classes: Uint32Array;
}

/**
 * One grid over all instances of the members, not one a member (a member's grid of instances spread over the scene is
 * about a cell an instance: the yeast cell's 31.5 M instances, 24.7 M cells split from its 5,639 members' grids). The
 * box the members' `bounds` (each holding its instances' spheres), the cells' side so that a cell holds
 * `cellInstances` instances on average over the box. A member's instances go into the grid of its size class, its
 * cells at least the member's sphere's diameter (`radii`: its invariant sphere's radius by its largest scale), so one
 * large member does not widen the cells of the small ones. Members without instances (`counts`) count for neither
 */
export function getSceneGrid(
  bounds: readonly Sphere3D[],
  radii: ArrayLike<number>,
  counts: ArrayLike<number>,
  instanceCount: number,
  cellInstances: number,
): SceneGrid {
  const segmentCount = bounds.length;
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  for (let s = 0; s < segmentCount; ++s) {
    if (counts[s] === 0) continue;
    const { center, radius } = bounds[s];
    for (let a = 0; a < 3; ++a) {
      lo[a] = Math.min(lo[a], center[a] - radius);
      hi[a] = Math.max(hi[a], center[a] + radius);
    }
  }
  // no instances: an empty box at the origin
  if (lo[0] > hi[0]) {
    lo.fill(0);
    hi.fill(0);
  }
  const ex = Math.max(hi[0] - lo[0], MinSceneCell);
  const ey = Math.max(hi[1] - lo[1], MinSceneCell);
  const ez = Math.max(hi[2] - lo[2], MinSceneCell);
  const side = Math.max(MinSceneCell, Math.cbrt((ex * ey * ez) / Math.max(1, instanceCount / cellInstances)));

  const classes = new Uint32Array(segmentCount);
  const used = new Array<boolean>(SceneCellClasses).fill(false);
  for (let s = 0; s < segmentCount; ++s) {
    const d = 2 * radii[s];
    const k = d <= side ? 0 : Math.min(SceneCellClasses - 1, Math.ceil(Math.log2(d / side)));
    classes[s] = k;
    if (counts[s] > 0) used[k] = true;
  }
  const grids = Mat4.zero();
  const firsts = Vec4();
  let cellCount = 0;
  for (let k = 0; k < SceneCellClasses; ++k) {
    firsts[k] = cellCount;
    if (!used[k]) continue;
    let l = side * 2 ** k;
    const product = Math.ceil(ex / l) * Math.ceil(ey / l) * Math.ceil(ez / l);
    if (product > MaxSceneCells) l *= Math.cbrt(product / MaxSceneCells) * 1.01;
    const nx = Math.max(1, Math.ceil(ex / l));
    const ny = Math.max(1, Math.ceil(ey / l));
    const nz = Math.max(1, Math.ceil(ez / l));
    grids[k * 4] = 1 / l;
    grids[k * 4 + 1] = nx;
    grids[k * 4 + 2] = ny;
    grids[k * 4 + 3] = nz;
    cellCount += nx * ny * nz;
  }
  return { min: Vec4.create(lo[0], lo[1], lo[2], 0), grids, firsts, cellCount, classes };
}

/** What a build reads and writes: the culling's buffers */
export interface SceneCellTargets {
  readonly instanceCount: number;
  readonly segmentData: ValueCell<StorageBuffer>;
  readonly data: ValueCell<StorageBuffer>;
  readonly instanceCells: ValueCell<StorageBuffer>;
  readonly cellSpheres: ValueCell<StorageBuffer>;
}

/** Builds scene cells on the GPU (`SceneCellBinKernel`), its kernels prepared (`SceneCellKernels`) */
export interface SceneCellBuilder {
  /** the grid and buffers of the next build: `instanceCells` and `cellSpheres` large enough for them */
  set: (grid: SceneGrid, targets: SceneCellTargets) => void;
  /** the boxes zeroed: outside a pass, before `dispatch` */
  clear: () => void;
  /** the build, in an open compute pass; `dispatchItems` a kernel over a count of items (in rows past the limit) */
  dispatch: (dispatchItems: (kernel: ComputeKernel, n: number) => void) => void;
  destroy: () => void;
}

export function createSceneCellBuilder(compute: ComputeContext): SceneCellBuilder {
  const buffer = (label: string) => ValueCell.create(compute.createBuffer(16, label));
  const cellClasses = buffer('scene cell classes');
  const cellBoxes = buffer('scene cell boxes');
  const ensure = (cell: ValueCell<StorageBuffer>, bytes: number, label: string) => {
    if (cell.ref.value.byteLength >= bytes) return;
    cell.ref.value.destroy();
    ValueCell.update(cell, compute.createBuffer(bytes, label));
  };
  const binValues = {
    uInstanceCount: ValueCell.create(0),
    uSceneMin: ValueCell.create(Vec4()),
    uSceneGrids: ValueCell.create(Mat4.zero()),
    uSceneFirsts: ValueCell.create(Vec4()),
    segmentData: ValueCell.create(compute.createBuffer(16)),
    cellClasses,
    data: ValueCell.create(compute.createBuffer(16)),
    instanceCells: ValueCell.create(compute.createBuffer(16)),
    cellBoxes,
  };
  // placeholders until `set` binds the culling's buffers
  const placeholders = [binValues.segmentData.ref.value, binValues.data.ref.value, binValues.instanceCells.ref.value];
  const spheresValues = {
    uCellCount: ValueCell.create(0),
    cellBoxes,
    cellSpheres: ValueCell.create(compute.createBuffer(16)),
  };
  placeholders.push(spheresValues.cellSpheres.ref.value);
  const bin = compute.createKernel(SceneCellBinKernel, SceneCellBinSchema, binValues);
  const spheres = compute.createKernel(SceneCellSpheresKernel, SceneCellSpheresSchema, spheresValues);
  let instanceCount = 0;
  let cellCount = 0;
  const bind = (cell: ValueCell<StorageBuffer>, from: ValueCell<StorageBuffer>) =>
    ValueCell.updateIfChanged(cell, from.ref.value);

  return {
    set: (grid, targets) => {
      instanceCount = targets.instanceCount;
      cellCount = grid.cellCount;
      ensure(cellClasses, grid.classes.byteLength, 'scene cell classes');
      cellClasses.ref.value.write(grid.classes);
      ensure(cellBoxes, cellCount * SceneCellBoxWords * 4, 'scene cell boxes');
      ValueCell.updateIfChanged(binValues.uInstanceCount, instanceCount);
      ValueCell.update(binValues.uSceneMin, Vec4.copy(binValues.uSceneMin.ref.value, grid.min));
      ValueCell.update(binValues.uSceneGrids, Mat4.copy(binValues.uSceneGrids.ref.value, grid.grids));
      ValueCell.update(binValues.uSceneFirsts, Vec4.copy(binValues.uSceneFirsts.ref.value, grid.firsts));
      ValueCell.updateIfChanged(spheresValues.uCellCount, cellCount);
      bind(binValues.segmentData, targets.segmentData);
      bind(binValues.data, targets.data);
      bind(binValues.instanceCells, targets.instanceCells);
      bind(spheresValues.cellSpheres, targets.cellSpheres);
    },
    clear: () => cellBoxes.ref.value.clear(0, cellCount * SceneCellBoxWords * 4),
    dispatch: (dispatchItems) => {
      if (instanceCount > 0) dispatchItems(bin, instanceCount);
      if (cellCount > 0) dispatchItems(spheres, cellCount);
    },
    destroy: () => {
      bin.destroy();
      spheres.destroy();
      cellClasses.ref.value.destroy();
      cellBoxes.ref.value.destroy();
      for (const b of placeholders) b.destroy();
    },
  };
}
