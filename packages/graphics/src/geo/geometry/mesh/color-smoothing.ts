/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { setColorType } from '../color-data.js';
import { setEmissiveType } from '../emissive-data.js';
import { setOverpaintType } from '../overpaint-data.js';
import { setSubstanceType } from '../substance-data.js';
import { setTransparencyType } from '../transparency-data.js';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';
import {
  createDataArray,
  type DataArray,
  type GpuDataWriter,
  getDataTextureImage,
  getSharedZeros,
  type TextureImage,
} from '@molstar/graphics/gl/renderable/util';
import { computeLayerSampling, isComputeLayerSamplingSupported } from '@molstar/graphics/gl/compute/layer-sampling';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Texture } from '@molstar/graphics/gl/texture';
import { Box3D, type Sphere3D } from '@molstar/core/math/geometry';
import { lerp } from '@molstar/core/math/interpolate';
import { Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { getVolumeTexture2dLayout } from '@molstar/graphics/geo/util';
import { ValueCell } from '@molstar/core/util';
import { getInstanceMatrices } from '../transform-data.js';
import { type GpuMeshVertices, getGpuMeshVertices } from './gpu-vertices.js';
import {
  computeMeshColorSmoothing,
  isComputeMeshColorSmoothingSupported,
} from '@molstar/graphics/gl/compute/color-smoothing';
import { createScratch } from '@molstar/graphics/gl/scratch';

/** The colors of the GPU splat (`computeMeshColorSmoothing`) by format, a texture each per context */
const ScratchTextures = createScratch<Texture>('mesh/color-smoothing colors');
/** The grids of the GPU splat: rgba for every layer */
const RgbaGrids = new WeakSet<Texture>();

interface ColorSmoothingInput {
  vertexCount: number;
  instanceCount: number;
  groupCount: number;
  transformBuffer: Float32Array;
  instanceBuffer: Uint32Array;
  positionBuffer: Float32Array;
  groupBuffer: Uint32Array;
  colorData: DataArray<Uint8Array>;
  colorType: 'group' | 'groupInstance';
  boundingSphere: Sphere3D;
  invariantBoundingSphere: Sphere3D;
  itemSize: 4 | 1;
  /**
   * the mesh's vertices on the GPU (`getGpuMeshVertices`): with a context of compute kernels for meshes (WebGPU), the
   * splat on the GPU (`computeMeshColorSmoothing`) into the same grid, the CPU's otherwise
   */
  gpuVertices?: GpuMeshVertices;
}

export type ColorSmoothingOptions = {
  resolution: number;
  stride: number;
};

/**
 * Grow-only scratch for the accumulators of `calcMeshColorSmoothing` (`data` with up to 4 values per cell, `count`),
 * all-zero between calls: the normalization clears what the splat wrote. The calculation is synchronous, so one set
 * serves every call; grids above `MaxScratchCells` (20 MB) get their own arrays so a rare large grid is not retained.
 */
const MaxScratchCells = 1 << 20;
const Scratch = { data: new Float32Array(0), count: new Float32Array(0), dirty: false };

function getAccumulators(cellCount: number, itemSize: number) {
  if (cellCount > MaxScratchCells)
    return { data: new Float32Array(cellCount * itemSize), count: new Float32Array(cellCount) };
  if (Scratch.count.length < cellCount) {
    const n = Math.min(MaxScratchCells, Math.max(cellCount, 2 * Scratch.count.length));
    Scratch.data = new Float32Array(n * 4);
    Scratch.count = new Float32Array(n);
  } else if (Scratch.dirty) {
    // a previous call threw before clearing
    Scratch.data.fill(0);
    Scratch.count.fill(0);
  }
  Scratch.dirty = true;
  return { data: Scratch.data, count: Scratch.count };
}

function releaseAccumulators(data: Float32Array) {
  if (data === Scratch.data) Scratch.dirty = false;
}

function getTextureSide(gridDim: Vec3) {
  const { width, height } = getVolumeTexture2dLayout(gridDim);
  return Math.max(width, height);
}

/**
 * Grid of the smoothing: `box` expanded by `1 + resolution`, two extra cells per axis, on the cell lattice of
 * `latticeBox` (the bounding-sphere box) expanded alike. A cell then has the same position and value, up to rounding,
 * as on the grid of `latticeBox`: its value depends only on the vertices within 2 cells of it, which are splatted
 * wherever the grid ends, and the vertices, which are what is interpolated, lie in `box` (as far as `latticeBox` holds
 * them). With `box` equal to `latticeBox`, this is the grid of `latticeBox` (as `getGridDim` in
 * `texture-mesh/color-smoothing.ts`). Vertex animation (wiggle, tumble) gets no slack beyond the padding, as on the
 * faces where the vertices touch `latticeBox`.
 */
function getGrid(latticeBox: Box3D, box: Box3D, resolution: number) {
  const pad = 1 + resolution;
  const scaleFactor = 1 / resolution;
  const min = Vec3();
  const gridDim = Vec3();
  for (let k = 0; k < 3; ++k) {
    const latticeMin = latticeBox.min[k] - pad;
    min[k] = latticeMin + Math.floor((box.min[k] - latticeBox.min[k]) * scaleFactor) * resolution;
    gridDim[k] = Math.ceil((box.max[k] + pad) * scaleFactor - min[k] * scaleFactor) + 2;
  }
  return { min, gridDim };
}

const tmpTransform = Mat4();
const tmpBox = Box3D();

/**
 * Box of the vertices (of all instances for instance types, from the transformed corners of the untransformed box
 * `vertexBox`), within `bound`; `bound` when there are no vertices
 */
function getVertexBox(
  input: ColorSmoothingInput,
  vertexBox: Box3D | undefined,
  instanceCount: number,
  isInstanceType: boolean,
  bound: Box3D,
) {
  if (!vertexBox) return Box3D.clone(bound);

  let box = Box3D.clone(vertexBox);
  if (isInstanceType) {
    const instances = Box3D.setEmpty(Box3D());
    for (let i = 0; i < instanceCount; ++i) {
      Mat4.fromArray(tmpTransform, input.transformBuffer, i * 16);
      Box3D.addBox3D(instances, Box3D.transform(tmpBox, box, tmpTransform));
    }
    box = instances;
  }

  // both boxes hold the vertices: take the smaller extent per axis (an empty intersection means `bound` misses
  // vertices, where the former grid cut them off too)
  Box3D.intersect(box, box, bound);
  return Box3D.isEmpty(box) ? Box3D.clone(bound) : box;
}

/** Splat radius, in cells */
const SplatRadius = 2;

/**
 * Cells each instance's vertices can splat into or are read from, 6 numbers per instance (min x, y, z and the
 * exclusive max x, y, z), not clamped to the grid: the transformed corners of `vertexBox` on the grid, the splat
 * radius beyond, and one more cell each side, since the corners and the vertices are transformed separately and may
 * round differently
 */
function getInstanceCells(
  transformBuffer: Float32Array,
  vertexBox: Box3D,
  instanceCount: number,
  isInstanceType: boolean,
  min: Vec3,
  scaleFactor: number,
) {
  const cells = new Float64Array(instanceCount * 6);
  for (let i = 0; i < instanceCount; ++i) {
    let box = vertexBox;
    if (isInstanceType) {
      Mat4.fromArray(tmpTransform, transformBuffer, i * 16);
      box = Box3D.transform(tmpBox, vertexBox, tmpTransform);
    }
    for (let k = 0; k < 3; ++k) {
      cells[i * 6 + k] = Math.floor((box.min[k] - min[k]) * scaleFactor) - SplatRadius - 1;
      cells[i * 6 + 3 + k] = Math.floor((box.max[k] - min[k]) * scaleFactor) + SplatRadius + 3;
    }
  }
  return cells;
}

/** Whether a vertex of a group in `groupBuffer` reads past the colors of its instance */
function hasGroupOutOfRange(groupBuffer: Uint32Array, vertexCount: number, groupCount: number) {
  for (let j = 0; j < vertexCount; ++j) {
    const g = groupBuffer[j];
    if (!(g >= 0 && g < groupCount)) return true;
  }
  return false;
}

/**
 * Whether any value of the bytes from `start` to `end` of `colors` is non-zero, read 4 bytes at a time through `words`
 * (a view of `colors`) when given; a range (partly) beyond `colors` counts as non-zero
 */
function isRangeReached(colors: Uint8Array, words: Uint32Array | undefined, start: number, end: number) {
  if (!(start >= 0 && end <= colors.length)) return true;
  let k = start;
  if (words) {
    for (; k < end && k % 4 !== 0; ++k) if (colors[k] !== 0) return true;
    const wordEnd = end >> 2;
    for (let w = k >> 2; w < wordEnd; ++w) if (words[w] !== 0) return true;
    k = Math.max(k, wordEnd * 4);
  }
  for (; k < end; ++k) if (colors[k] !== 0) return true;
  return false;
}

/**
 * Whether instance row `row` of the layer `colors` is reached (any value non-zero), with the `itemsBefore` items before
 * the row, which vertices of negative groups read; a row (partly) beyond `colors` counts as reached
 */
export function getRowReachedTest(colors: Uint8Array, groupCount: number, itemSize: number, itemsBefore = 0) {
  const words =
    colors.byteOffset % 4 === 0 ? new Uint32Array(colors.buffer, colors.byteOffset, colors.length >> 2) : undefined;
  return (row: number) => {
    const start = row * groupCount * itemSize;
    return isRangeReached(colors, words, start - itemsBefore * itemSize, start + groupCount * itemSize);
  };
}

/** Whether the cells of instance `i` overlap the box at `o` of `boxes` (6 numbers per box, as `cells`) */
function overlaps(cells: Float64Array, i: number, boxes: Float64Array, o: number) {
  const c = i * 6;
  return (
    cells[c] < boxes[o + 3] &&
    cells[c + 3] > boxes[o] &&
    cells[c + 1] < boxes[o + 4] &&
    cells[c + 4] > boxes[o + 1] &&
    cells[c + 2] < boxes[o + 5] &&
    cells[c + 5] > boxes[o + 2]
  );
}

/** Above this many reached instances, the others are tested against the region only, not each reached instance */
const MaxReachedTests = 256;

/**
 * The instances to splat (perf-audit 13/P1, `.v6/perf-audit/13-large-scale.md`): an instance whose layer is non-zero
 * ("reached", `isReached`), and an instance whose cells overlap the cells of a reached one. The region, where the
 * accumulators are, is the box of the reached instances' cells, clamped to the grid. `vertexBox` holds the untransformed
 * vertices; the grid is at `min`, `scaleFactor` cells per unit, `gridDim` cells.
 *
 * Exact: only cells of reached instances can become non-zero, and every vertex that adds to such a cell is splatted
 * in the former order, so these cells are bit-identical to splatting every instance; every other cell gets zeros or
 * nothing, so it is 0, as before. A reached instance with non-finite cells (a degenerate transform) splats every
 * instance over the whole grid, as before. Also used by the GPU path of texture meshes (`calcTextureMeshColorSmoothing`).
 *
 * `read` (the CPU path only): the instances to interpolate, those splatted and those whose vertices are not within the
 * grid (they read cells of other slices, as before); the others read only cells of their own, which are 0.
 */
export function getSplatInstances(
  transformBuffer: Float32Array,
  vertexBox: Box3D,
  instanceCount: number,
  isInstanceType: boolean,
  isReached: (i: number) => boolean,
  gridDim: Vec3,
  min: Vec3,
  scaleFactor: number,
  withRead: boolean,
) {
  const cells = getInstanceCells(transformBuffer, vertexBox, instanceCount, isInstanceType, min, scaleFactor);

  const splat = new Uint8Array(instanceCount);
  const reached: number[] = [];
  let finite = true;
  for (let i = 0; i < instanceCount; ++i) {
    if (!isReached(i)) continue;
    splat[i] = 1;
    reached.push(i);
    for (let k = 0; k < 6; ++k) if (!Number.isFinite(cells[i * 6 + k])) finite = false;
  }

  const region = new Float64Array([0, 0, 0, gridDim[0], gridDim[1], gridDim[2]]);
  if (!finite) {
    splat.fill(1);
    return { region, splat, read: withRead ? splat : undefined };
  }
  if (reached.length === 0) return { region: undefined, splat, read: withRead ? splat : undefined };

  region.set([Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity]);
  for (const i of reached) {
    for (let k = 0; k < 3; ++k) {
      region[k] = Math.min(region[k], cells[i * 6 + k]);
      region[3 + k] = Math.max(region[3 + k], cells[i * 6 + 3 + k]);
    }
  }
  for (let k = 0; k < 3; ++k) {
    region[k] = Math.max(0, region[k]);
    region[3 + k] = Math.min(gridDim[k], region[3 + k]);
  }
  const isEmpty = !(region[0] < region[3] && region[1] < region[4] && region[2] < region[5]);
  if (isEmpty) return { region: undefined, splat: splat.fill(0), read: withRead ? splat : undefined };

  // the neighbours: the region first, then the reached instances (unless there are many, the region alone)
  const testEach = reached.length <= MaxReachedTests;
  for (let i = 0; i < instanceCount; ++i) {
    if (splat[i] || !overlaps(cells, i, region, 0)) continue;
    if (!testEach) {
      splat[i] = 1;
      continue;
    }
    for (const r of reached) {
      if (overlaps(cells, i, cells, r * 6)) {
        splat[i] = 1;
        break;
      }
    }
  }
  if (!withRead) return { region, splat, read: undefined };

  const read = new Uint8Array(instanceCount);
  for (let i = 0; i < instanceCount; ++i) {
    // the cells read (the trilinear neighbours) are the instance's cells without the splat radius
    const o = i * 6;
    const r = SplatRadius;
    const inGrid =
      cells[o] + r >= 0 &&
      cells[o + 1] + r >= 0 &&
      cells[o + 2] + r >= 0 &&
      cells[o + 3] - r <= gridDim[0] &&
      cells[o + 4] - r <= gridDim[1] &&
      cells[o + 5] - r <= gridDim[2];
    read[i] = splat[i] || !inGrid ? 1 : 0;
  }
  return { region, splat, read };
}

interface SplatContext {
  positionBuffer: Float32Array;
  transformBuffer: Float32Array;
  groupBuffer: Uint32Array;
  colors: Uint8Array;
  isInstanceType: boolean;
  groupCount: number;
  itemSize: number;
  min: Vec3;
  scaleFactor: number;
  /** the cells of the accumulators (min x, y, z, exclusive max x, y, z), x fastest */
  region: Float64Array;
  data: Float32Array;
  count: Float32Array;
  v: Vec3;
}

/**
 * Splats vertex `j` of instance `i` (colors of row `instanceIndex`) into the accumulators within the region. A
 * function of its own, called per vertex, so that V8 optimizes it early: as part of the caller, which runs once per
 * smoothing, the loops ran on on-stack-replaced code (a cold one-chain overpaint on 1tup ×864 at `high`: 42 → 31 ms,
 * 50 → 25 MB of JS heap allocated)
 */
function splatVertex(ctx: SplatContext, i: number, instanceIndex: number, j: number) {
  const { colors, itemSize, region, data, count, v } = ctx;
  Vec3.fromArray(v, ctx.positionBuffer, j * 3);
  if (ctx.isInstanceType) Vec3.transformMat4Offset(v, v, ctx.transformBuffer, 0, 0, i * 16);
  Vec3.sub(v, v, ctx.min);
  Vec3.scale(v, v, ctx.scaleFactor);
  const vx = v[0];
  const vy = v[1];
  const vz = v[2];

  // vertex mapped to grid
  const x = Math.floor(vx);
  const y = Math.floor(vy);
  const z = Math.floor(vz);

  const p = SplatRadius;
  const p2 = p * p;
  const rx0 = region[0];
  const ry0 = region[1];
  const rz0 = region[2];
  const rx = region[3] - rx0;
  const ry = region[4] - ry0;

  // Extents of grid to consider for this atom, within the region
  const begX = Math.max(rx0, x - p);
  const begY = Math.max(ry0, y - p);
  const begZ = Math.max(rz0, z - p);

  // Add two to these points:
  // - x, y, z are floor'd values so this ensures coverage
  // - these are loop limits (exclusive)
  const endX = Math.min(region[3], x + p + 2);
  const endY = Math.min(region[4], y + p + 2);
  const endZ = Math.min(region[5], z + p + 2);
  if (begX >= endX || begY >= endY || begZ >= endZ) return;

  // group colors
  const ci = (instanceIndex * ctx.groupCount + ctx.groupBuffer[j]) * itemSize;
  const c0 = colors[ci];
  const c1 = itemSize === 4 ? colors[ci + 1] : 0;
  const c2 = itemSize === 4 ? colors[ci + 2] : 0;
  const c3 = itemSize === 4 ? colors[ci + 3] : 0;

  // z/y/x with x innermost; each cell gets at most one contribution per vertex, so the per-cell summation order is
  // unchanged. `dz2 > p2` and `dy2 + dz2 > p2` only skip cells whose distance below exceeds p as well (rounded addition
  // of non-negative terms is monotone), and `d2` keeps the original association so `d` is bit-identical; cells at
  // `s == 0` (`d2 > p2` or `d == p`) would add 0, which changes nothing for finite colors. Skipping them also keeps the
  // scratch clean: every write makes `count` > 0, so the normalization clears it, even for a non-finite color (an
  // out-of-range group or instance).
  for (let zi = begZ; zi < endZ; ++zi) {
    const dz = zi - vz;
    const dz2 = dz * dz;
    if (dz2 > p2) continue;
    const zo = (zi - rz0) * ry;
    for (let yi = begY; yi < endY; ++yi) {
      const dy = yi - vy;
      const dy2 = dy * dy;
      if (dy2 + dz2 > p2) continue;
      const yo = (zo + yi - ry0) * rx - rx0;
      for (let xi = begX; xi < endX; ++xi) {
        const dx = xi - vx;
        const d2 = dx * dx + dy2 + dz2;
        if (d2 > p2) continue;

        const s = p - Math.sqrt(d2);
        if (s === 0) continue;
        const c = yo + xi;
        if (itemSize === 4) {
          const index = c * 4;
          data[index] += c0 * s;
          data[index + 1] += c1 * s;
          data[index + 2] += c2 * s;
          data[index + 3] += c3 * s;
        } else {
          data[c] += c0 * s;
        }
        count[c] += s;
      }
    }
  }
}

export function calcMeshColorSmoothing(
  input: ColorSmoothingInput,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  texture?: Texture,
) {
  // a texture of another context (e.g. before a backend switch) is not reused
  if (texture && graphics && !graphics.resources.ownsTexture(texture)) texture = undefined;
  const { colorType, vertexCount, groupCount, positionBuffer, instanceBuffer, transformBuffer, groupBuffer, itemSize } =
    input;
  const { stride } = options;

  const isInstanceType = colorType.endsWith('Instance');
  const instanceCount = isInstanceType ? input.instanceCount : 1;
  const sphereBox = Box3D.fromSphere3D(Box3D(), isInstanceType ? input.boundingSphere : input.invariantBoundingSphere);
  const vertexBox = vertexCount > 0 ? Box3D.fromPositionArray(Box3D(), positionBuffer, vertexCount) : undefined;
  const box = getVertexBox(input, vertexBox, instanceCount, isInstanceType, sphereBox);
  // the grid is uploaded into a texture with a graphics context: coarsen as `getFittingResolution`, on this grid
  let resolution = options.resolution;
  if (graphics) {
    const { maxTextureSize } = graphics.resources.limits;
    while (getTextureSide(getGrid(sphereBox, box, resolution).gridDim) > maxTextureSize) resolution *= 1.1;
  }
  // the grid spans every instance, also where the layer is zero: the shaders read it without a bounds check
  // (`texture3dFrom2dLinear` in WebGL and WebGPU), so a point beyond a smaller grid would read another slice
  const { min, gridDim } = getGrid(sphereBox, box, resolution);
  const scaleFactor = 1 / resolution;

  const { width, height } = getVolumeTexture2dLayout(gridDim);

  const grid = new Uint8Array(width * height * itemSize);
  const textureImage: TextureImage<Uint8Array> = { array: grid, width, height, filter: 'linear' };

  const colors = input.colorData.array;

  const zOff = getSliceOffsets(gridDim, width);

  // - use reordered index for access from GPU
  // - use serial index for access from CPU
  const instanceRow = (i: number) => (graphics && isInstanceType ? instanceBuffer[i] : i);
  const allReached = hasGroupOutOfRange(groupBuffer, vertexCount, groupCount);
  const isRowReached = getRowReachedTest(colors, groupCount, itemSize);
  const instances = vertexBox
    ? getSplatInstances(
        transformBuffer,
        vertexBox,
        instanceCount,
        isInstanceType,
        (i) => allReached || isRowReached(instanceRow(i)),
        gridDim,
        min,
        scaleFactor,
        !graphics,
      )
    : undefined;

  // the splat on the GPU: the same grid, instances and region, the sums of the kernels (`computeMeshColorSmoothing`)
  const gpuVertices = input.gpuVertices;
  if (
    graphics &&
    gpuVertices &&
    gpuVertices.compute === graphics.compute &&
    isComputeMeshColorSmoothingSupported(graphics)
  ) {
    // an rgba texture also for single-channel layers: `r8unorm` is no storage format in core WebGPU (the readers read
    // red); a single-channel one of the CPU's (`r`) not reused
    if (texture && itemSize === 1 && !RgbaGrids.has(texture)) texture = undefined;
    if (!texture) {
      texture = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
      RgbaGrids.add(texture);
    }
    const colorTexture = ScratchTextures.getOrCreate(
      graphics,
      () => graphics.resources.texture('image-uint8', itemSize === 4 ? 'rgba' : 'r', 'ubyte', 'nearest'),
      itemSize === 4 ? 'rgba' : 'r',
    );
    colorTexture.load(getDataTextureImage(input.colorData, itemSize));
    let splatCount = 0;
    const splat = instances?.splat;
    for (let i = 0; i < instanceCount; ++i) if (!splat || splat[i]) splatCount += 1;
    const transforms = new Float32Array(Math.max(1, splatCount) * 16);
    const rows = new Uint32Array(Math.max(1, splatCount));
    for (let i = 0, k = 0; i < instanceCount; ++i) {
      if (splat && !splat[i]) continue;
      if (isInstanceType) transforms.set(transformBuffer.subarray(i * 16, i * 16 + 16), k * 16);
      else Mat4.toArray(Mat4.id, transforms, k * 16);
      rows[k] = instanceRow(i);
      k += 1;
    }
    const empty = [0, 0, 0, 0, 0, 0];
    computeMeshColorSmoothing(
      graphics.compute,
      {
        vertexCount,
        stride,
        sampleCount: Math.ceil(vertexCount / stride),
        instanceCount: splatCount,
        groupCount,
        transformBuffer: transforms,
        instanceBuffer: rows,
        colorData: colorTexture,
        colorType,
        positions: gpuVertices.positions,
        groups: gpuVertices.groups,
        vertexOffset: 0,
      },
      {
        gridDim,
        texCols: width / gridDim[0],
        width,
        height,
        bboxMin: min,
        bboxSize: Vec3.scale(Vec3(), gridDim, resolution),
        resolution,
        region: vertexBox && instances?.region ? instances.region : empty,
      },
      texture,
    );
    return {
      kind: 'volume' as const,
      texture,
      gridTexDim: Vec2.create(width, height),
      gridDim,
      gridTransform: Vec4.create(min[0], min[1], min[2], scaleFactor),
      type: isInstanceType ? ('volumeInstance' as const) : ('volume' as const),
      splatCount,
      reached: instances?.read,
    };
  }

  /** the instances splatted (perf-audit 13/P1) */
  let splatCount = 0;
  if (instances?.region) {
    const { region, splat } = instances;
    // the accumulators span the region only, x fastest
    const [rx0, ry0, rz0] = region;
    const rx = region[3] - rx0;
    const ry = region[4] - ry0;
    const rz = region[5] - rz0;
    const cellCount = rx * ry * rz;
    const { data, count } = getAccumulators(cellCount, itemSize);

    const ctx: SplatContext = {
      positionBuffer,
      transformBuffer,
      groupBuffer,
      colors,
      isInstanceType,
      groupCount,
      itemSize,
      min,
      scaleFactor,
      region,
      data,
      count,
      v: Vec3(),
    };
    for (let i = 0; i < instanceCount; ++i) {
      if (!splat[i]) continue;
      splatCount += 1;
      const instanceIndex = instanceRow(i);
      for (let j = 0; j < vertexCount; j += stride) splatVertex(ctx, i, instanceIndex, j);
    }

    const rxEnd = region[3];
    const ryEnd = region[4];
    const rzEnd = region[5];
    // cells without a contribution stay 0 in `grid` (formerly `Math.round(0 / 0)`, also 0 in a Uint8Array); a cell
    // with `data` > 0 has `count` > 0, so clearing the counted cells leaves the accumulators all-zero for the next call
    for (let zi = rz0; zi < rzEnd; ++zi) {
      for (let yi = ry0; yi < ryEnd; ++yi) {
        const l0 = ((zi - rz0) * ry + yi - ry0) * rx - rx0;
        const g0 = zOff[zi] + yi * width;
        for (let xi = rx0; xi < rxEnd; ++xi) {
          const l = l0 + xi;
          const c = count[l];
          if (c === 0) continue;
          const g = g0 + xi;
          if (itemSize === 4) {
            const ls = l * 4;
            const gs = g * 4;
            grid[gs] = Math.round(data[ls] / c);
            grid[gs + 1] = Math.round(data[ls + 1] / c);
            grid[gs + 2] = Math.round(data[ls + 2] / c);
            grid[gs + 3] = Math.round(data[ls + 3] / c);
            data[ls] = 0;
            data[ls + 1] = 0;
            data[ls + 2] = 0;
            data[ls + 3] = 0;
          } else {
            grid[g] = Math.round(data[l] / c);
            data[l] = 0;
          }
          count[l] = 0;
        }
      }
    }
    releaseAccumulators(data);
  }

  const gridTexDim = Vec2.create(width, height);
  const gridTransform = Vec4.create(min[0], min[1], min[2], scaleFactor);
  const type = isInstanceType ? ('volumeInstance' as const) : ('volume' as const);

  if (graphics) {
    if (!texture) {
      const format = itemSize === 4 ? 'rgba' : 'r';
      texture = graphics.resources.texture('image-uint8', format, 'ubyte', 'linear');
    }
    texture.load(textureImage);

    // the instances reached (`read`, of instance types), for samplers of the grid per vertex on the GPU
    return {
      kind: 'volume' as const,
      texture,
      gridTexDim,
      gridDim,
      gridTransform,
      type,
      splatCount,
      reached: instances?.read,
    };
  } else {
    const interpolated = getTrilinearlyInterpolated({
      vertexCount,
      instanceCount,
      transformBuffer,
      positionBuffer,
      colorType: type,
      grid,
      gridDim,
      gridTexDim,
      gridTransform,
      vertexStride: 3,
      colorStride: itemSize,
      outputStride: itemSize,
      instances: instances?.read,
    });

    return {
      kind: 'vertex' as const,
      splatCount,
      texture: interpolated,
      type: isInstanceType ? 'vertexInstance' : 'vertex',
    };
  }
}

/**
 * Offset of slice `z` in the 2d layout of a grid (`getVolumeTexture2dLayout`), in cells: cell (x, y, z) is at
 * `offsets[z] + y * width + x`
 */
function getSliceOffsets(gridDim: Vec3, width: number) {
  const [xn, yn, zn] = gridDim;
  const offsets = new Float64Array(zn);
  for (let z = 0; z < zn; ++z) {
    const column = Math.floor(((z * xn) % width) / xn);
    const row = Math.floor((z * xn) / width);
    offsets[z] = row * yn * width + column * xn;
  }
  return offsets;
}

//

interface ColorInterpolationInput {
  vertexCount: number;
  instanceCount: number;
  transformBuffer: Float32Array;
  positionBuffer: Float32Array;
  colorType: 'volumeInstance' | 'volume';
  grid: Uint8Array; // 2d layout
  gridTexDim: Vec2;
  gridDim: Vec3;
  gridTransform: Vec4;
  vertexStride: 3 | 4;
  colorStride: 1 | 3 | 4;
  outputStride: 1 | 3 | 4;
  itemOffset?: 0 | 1 | 2 | 3;
  /** the instances to interpolate (non-zero), the others stay 0; all when not given */
  instances?: Uint8Array;
}

export function getTrilinearlyInterpolated(input: ColorInterpolationInput): DataArray<Uint8Array> {
  const {
    vertexCount,
    positionBuffer,
    transformBuffer,
    grid,
    gridDim,
    gridTexDim,
    gridTransform,
    vertexStride,
    colorStride,
  } = input;

  const itemOffset = input.itemOffset || 0;
  const outputStride = input.outputStride;
  if (outputStride + itemOffset > colorStride) {
    throw new Error('outputStride + itemOffset must NOT be larger than colorStride');
  }

  const isInstanceType = input.colorType.endsWith('Instance');
  const instanceCount = isInstanceType ? input.instanceCount : 1;
  const image = createDataArray(Math.max(1, instanceCount * vertexCount), outputStride, Uint8Array);
  const { array } = image;

  const [xn, yn, zn] = gridDim;
  const width = gridTexDim[0];
  const minX = gridTransform[0];
  const minY = gridTransform[1];
  const minZ = gridTransform[2];
  const scaleFactor = gridTransform[3];

  const zOff = getSliceOffsets(gridDim, width);
  /** a slice's offset in the 2d layout: vertices outside the grid (not expected) keep their former, layout-derived one */
  function sliceOffset(z: number) {
    return z >= 0 && z < zn
      ? zOff[z]
      : Math.floor(((z * xn) % width) / xn) * xn + Math.floor((z * xn) / width) * yn * width;
  }

  // scalars, not `Vec3`s, in this loop over every vertex of every instance (array destructuring and calls per vertex:
  // a smoothed layer set on all 420 units of 60 copies of 7qpd 3.67 to 2.75 s a frame, WebGPU, its splatting the rest); the
  // same double arithmetic in the same order as `Vec3.transformMat4Offset`, `sub` and `scale`, so the same bytes
  const { instances } = input;
  for (let i = 0; i < instanceCount; ++i) {
    if (instances && !instances[i]) continue;
    const t = i * 16;
    for (let j = 0; j < vertexCount; ++j) {
      const p = j * vertexStride;
      let x = positionBuffer[p];
      let y = positionBuffer[p + 1];
      let z = positionBuffer[p + 2];
      if (isInstanceType) {
        const m = transformBuffer;
        const w = 1 / (m[3 + t] * x + m[7 + t] * y + m[11 + t] * z + m[15 + t] || 1.0);
        const tx = (m[0 + t] * x + m[4 + t] * y + m[8 + t] * z + m[12 + t]) * w;
        const ty = (m[1 + t] * x + m[5 + t] * y + m[9 + t] * z + m[13 + t]) * w;
        const tz = (m[2 + t] * x + m[6 + t] * y + m[10 + t] * z + m[14 + t]) * w;
        x = tx;
        y = ty;
        z = tz;
      }
      x = (x - minX) * scaleFactor;
      y = (y - minY) * scaleFactor;
      z = (z - minZ) * scaleFactor;

      const x0 = Math.floor(x);
      const y0 = Math.floor(y);
      const z0 = Math.floor(z);
      const x1 = Math.ceil(x);
      const y1 = Math.ceil(y);
      const z1 = Math.ceil(z);

      // the cell is 1 wide, or 0 on a grid plane, where the weight is 0: not 0/0 (perf-audit-2 08/B3)
      const xd = x - x0;
      const yd = y - y0;
      const zd = z - z0;

      const zo0 = sliceOffset(z0);
      const zo1 = sliceOffset(z1);
      const r00 = y0 * width;
      const r10 = y1 * width;
      const i000 = colorStride * (zo0 + r00 + x0) + itemOffset;
      const i100 = colorStride * (zo0 + r00 + x1) + itemOffset;
      const i001 = colorStride * (zo1 + r00 + x0) + itemOffset;
      const i101 = colorStride * (zo1 + r00 + x1) + itemOffset;
      const i010 = colorStride * (zo0 + r10 + x0) + itemOffset;
      const i110 = colorStride * (zo0 + r10 + x1) + itemOffset;
      const i011 = colorStride * (zo1 + r10 + x0) + itemOffset;
      const i111 = colorStride * (zo1 + r10 + x1) + itemOffset;

      const o = (i * vertexCount + j) * outputStride;

      for (let k = 0; k < outputStride; ++k) {
        const s000 = grid[i000 + k];
        const s100 = grid[i100 + k];
        const s001 = grid[i001 + k];
        const s101 = grid[i101 + k];
        const s010 = grid[i010 + k];
        const s110 = grid[i110 + k];
        const s011 = grid[i011 + k];
        const s111 = grid[i111 + k];

        const s00 = lerp(s000, s100, xd);
        const s01 = lerp(s001, s101, xd);
        const s10 = lerp(s010, s110, xd);
        const s11 = lerp(s011, s111, xd);

        const s0 = lerp(s00, s10, yd);
        const s1 = lerp(s01, s11, yd);

        array[o + k] = lerp(s0, s1, zd);
      }
    }
  }

  return image;
}

//

/**
 * A smoothed layer as `vertexInstance` data written on the GPU (`DataArray.gpu`): its grid (`smoothing`, splatted on the
 * CPU and uploaded) sampled per vertex and instance by `computeLayerSampling` into the storage its render item draws it
 * from; its array zeros (`getSharedZeros`). The positions those of the mesh on the GPU (`getGpuMeshVertices`, once a
 * version of them), resolved at each write; the transforms and reached flags the layer's.
 */
function createGpuLayerData(
  values: MeshValues,
  smoothing: { texture: Texture; gridDim: Vec3; gridTexDim: Vec2; gridTransform: Vec4; reached?: Uint8Array },
  itemSize: 1 | 4,
  graphics: GraphicsContext,
): DataArray<Uint8Array> {
  const compute = graphics.compute!;
  const vertexCount = values.uVertexCount.ref.value;
  const instanceCount = values.uInstanceCount.ref.value;
  const count = vertexCount * instanceCount;
  const matrices = getInstanceMatrices(values).subarray(0, instanceCount * 16);
  const transforms = compute.createBuffer(Math.max(64, matrices.byteLength), 'layer sampling transforms');
  transforms.write(matrices);
  const flags = new Uint32Array(Math.max(1, instanceCount));
  for (let i = 0; i < instanceCount; ++i) flags[i] = !smoothing.reached || smoothing.reached[i] ? 1 : 0;
  const reached = compute.createBuffer(flags.byteLength, 'layer sampling reached');
  reached.write(flags);
  const vertices = () => getGpuMeshVertices(values, compute).positions;
  const input = {
    vertexCount,
    instanceCount,
    positions: vertices(),
    vertexOffset: 0,
    transforms,
    instanceOffset: 0,
    reached,
    grid: smoothing.texture,
    gridDim: smoothing.gridDim,
    gridTexDim: smoothing.gridTexDim,
    gridTransform: smoothing.gridTransform,
    itemSize,
  };
  const byteLength = count * itemSize;
  const writer: GpuDataWriter = {
    compute,
    elementBytes: itemSize,
    write: (buffer, byteOffset) => {
      input.positions = vertices();
      computeLayerSampling(compute, input, buffer, byteOffset / itemSize);
    },
    read: async () => {
      // written into a buffer of its own: where the render item holds it may have moved since
      const out = compute.createBuffer(Math.max(4, Math.ceil(byteLength / 4) * 4), 'layer sampling read');
      input.positions = vertices();
      computeLayerSampling(compute, input, out, 0);
      const bytes = new Uint8Array(await out.readAsync(), 0, byteLength).slice();
      out.destroy();
      return bytes;
    },
    dispose: () => {
      transforms.destroy();
      reached.destroy();
      smoothing.texture.destroy();
    },
  };
  return { array: getSharedZeros(Uint8Array, byteLength), count, gpu: writer };
}

/**
 * The mesh's vertices on the GPU where `graphics` splats meshes there (WebGPU compute): the grid splatted by the kernels
 * (`computeMeshColorSmoothing`), the CPU's splat otherwise
 */
function gpuVerticesOf(values: MeshValues, graphics: GraphicsContext | undefined) {
  if (!graphics || !isComputeMeshColorSmoothingSupported(graphics)) return undefined;
  return getGpuMeshVertices(values, graphics.compute);
}

function isSupportedColorType(x: string): x is 'group' | 'groupInstance' {
  return x === 'group' || x === 'groupInstance';
}

export function applyMeshColorSmoothing(
  values: MeshValues,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedColorType(values.colorType.ref.value)) return;

  const smoothingData = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: getInstanceMatrices(values),
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: values.tColor.ref.value,
      colorType: values.colorType.ref.value,
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize: 4,
      gpuVertices: gpuVerticesOf(values, graphics),
    },
    options,
    graphics,
    colorTexture,
  );

  if (smoothingData.kind === 'volume') {
    setColorType(values, smoothingData.type);
    ValueCell.update(values.tColorGrid, smoothingData.texture);
    ValueCell.update(values.uColorGridTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uColorGridDim, smoothingData.gridDim);
    ValueCell.update(values.uColorGridTransform, smoothingData.gridTransform);
  } else if (smoothingData.kind === 'vertex') {
    setColorType(values, smoothingData.type);
    ValueCell.update(values.tColor, smoothingData.texture);
  }
}

function isSupportedOverpaintType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyMeshOverpaintSmoothing(
  values: MeshValues,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  colorTexture?: Texture,
  perVertexOnGpu?: boolean,
) {
  if (!isSupportedOverpaintType(values.overpaintType.ref.value)) return;

  // per vertex: sampled on the GPU where it computes (its grid uploaded), else on the CPU (no context, no grid)
  const gpu = !!perVertexOnGpu && !!graphics && isComputeLayerSamplingSupported(graphics);
  const smoothingData = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: getInstanceMatrices(values),
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: values.tOverpaint.ref.value,
      colorType: values.overpaintType.ref.value,
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize: 4,
      gpuVertices: gpuVerticesOf(values, gpu || !perVertexOnGpu ? graphics : undefined),
    },
    options,
    gpu || !perVertexOnGpu ? graphics : undefined,
    gpu ? undefined : colorTexture,
  );
  if (gpu && smoothingData.kind === 'volume') {
    setOverpaintType(values, 'vertexInstance');
    ValueCell.update(values.tOverpaint, createGpuLayerData(values, smoothingData, 4, graphics!));
    return;
  }
  if (smoothingData.kind === 'volume') {
    setOverpaintType(values, smoothingData.type);
    ValueCell.update(values.tOverpaintGrid, smoothingData.texture);
    ValueCell.update(values.uOverpaintGridTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uOverpaintGridDim, smoothingData.gridDim);
    ValueCell.update(values.uOverpaintGridTransform, smoothingData.gridTransform);
  } else if (smoothingData.kind === 'vertex') {
    setOverpaintType(values, smoothingData.type);
    ValueCell.update(values.tOverpaint, smoothingData.texture);
  }
}

function isSupportedTransparencyType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyMeshTransparencySmoothing(
  values: MeshValues,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  colorTexture?: Texture,
  perVertexOnGpu?: boolean,
) {
  if (!isSupportedTransparencyType(values.transparencyType.ref.value)) return;

  // per vertex: sampled on the GPU where it computes (its grid uploaded), else on the CPU (no context, no grid)
  const gpu = !!perVertexOnGpu && !!graphics && isComputeLayerSamplingSupported(graphics);
  const smoothingData = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: getInstanceMatrices(values),
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: values.tTransparency.ref.value,
      colorType: values.transparencyType.ref.value,
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize: 1,
      gpuVertices: gpuVerticesOf(values, gpu || !perVertexOnGpu ? graphics : undefined),
    },
    options,
    gpu || !perVertexOnGpu ? graphics : undefined,
    gpu ? undefined : colorTexture,
  );
  if (gpu && smoothingData.kind === 'volume') {
    setTransparencyType(values, 'vertexInstance');
    ValueCell.update(values.tTransparency, createGpuLayerData(values, smoothingData, 1, graphics!));
    return;
  }
  if (smoothingData.kind === 'volume') {
    setTransparencyType(values, smoothingData.type);
    ValueCell.update(values.tTransparencyGrid, smoothingData.texture);
    ValueCell.update(values.uTransparencyGridTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uTransparencyGridDim, smoothingData.gridDim);
    ValueCell.update(values.uTransparencyGridTransform, smoothingData.gridTransform);
  } else if (smoothingData.kind === 'vertex') {
    setTransparencyType(values, smoothingData.type);
    ValueCell.update(values.tTransparency, smoothingData.texture);
  }
}

function isSupportedEmissiveType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyMeshEmissiveSmoothing(
  values: MeshValues,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  colorTexture?: Texture,
  perVertexOnGpu?: boolean,
) {
  if (!isSupportedEmissiveType(values.emissiveType.ref.value)) return;

  // per vertex: sampled on the GPU where it computes (its grid uploaded), else on the CPU (no context, no grid)
  const gpu = !!perVertexOnGpu && !!graphics && isComputeLayerSamplingSupported(graphics);
  const smoothingData = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: getInstanceMatrices(values),
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: values.tEmissive.ref.value,
      colorType: values.emissiveType.ref.value,
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize: 1,
      gpuVertices: gpuVerticesOf(values, gpu || !perVertexOnGpu ? graphics : undefined),
    },
    options,
    gpu || !perVertexOnGpu ? graphics : undefined,
    gpu ? undefined : colorTexture,
  );
  if (gpu && smoothingData.kind === 'volume') {
    setEmissiveType(values, 'vertexInstance');
    ValueCell.update(values.tEmissive, createGpuLayerData(values, smoothingData, 1, graphics!));
    return;
  }
  if (smoothingData.kind === 'volume') {
    setEmissiveType(values, smoothingData.type);
    ValueCell.update(values.tEmissiveGrid, smoothingData.texture);
    ValueCell.update(values.uEmissiveGridTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uEmissiveGridDim, smoothingData.gridDim);
    ValueCell.update(values.uEmissiveGridTransform, smoothingData.gridTransform);
  } else if (smoothingData.kind === 'vertex') {
    setEmissiveType(values, smoothingData.type);
    ValueCell.update(values.tEmissive, smoothingData.texture);
  }
}

function isSupportedSubstanceType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyMeshSubstanceSmoothing(
  values: MeshValues,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  colorTexture?: Texture,
  perVertexOnGpu?: boolean,
) {
  if (!isSupportedSubstanceType(values.substanceType.ref.value)) return;

  // per vertex: sampled on the GPU where it computes (its grid uploaded), else on the CPU (no context, no grid)
  const gpu = !!perVertexOnGpu && !!graphics && isComputeLayerSamplingSupported(graphics);
  const smoothingData = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: getInstanceMatrices(values),
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: values.tSubstance.ref.value,
      colorType: values.substanceType.ref.value,
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize: 4,
      gpuVertices: gpuVerticesOf(values, gpu || !perVertexOnGpu ? graphics : undefined),
    },
    options,
    gpu || !perVertexOnGpu ? graphics : undefined,
    gpu ? undefined : colorTexture,
  );
  if (gpu && smoothingData.kind === 'volume') {
    setSubstanceType(values, 'vertexInstance');
    ValueCell.update(values.tSubstance, createGpuLayerData(values, smoothingData, 4, graphics!));
    return;
  }
  if (smoothingData.kind === 'volume') {
    setSubstanceType(values, smoothingData.type);
    ValueCell.update(values.tSubstanceGrid, smoothingData.texture);
    ValueCell.update(values.uSubstanceGridTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uSubstanceGridDim, smoothingData.gridDim);
    ValueCell.update(values.uSubstanceGridTransform, smoothingData.gridTransform);
  } else if (smoothingData.kind === 'vertex') {
    setSubstanceType(values, smoothingData.type);
    ValueCell.update(values.tSubstance, smoothingData.texture);
  }
}
