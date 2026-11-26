/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { getTrilinearlyInterpolated } from '@molstar/graphics/geo/geometry/mesh/color-smoothing';
import { computeLayerSampling, isComputeLayerSamplingSupported } from '@molstar/graphics/gl/compute/layer-sampling';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createWebGPUPlugin } from './webgpu-plugin.js';
import { WebGPU } from './webgpu-module.js';

const GridDim = Vec3.create(7, 5, 6);
/** slices in rows of three, as the smoothing's 2d atlas lays them out */
const Columns = 3;
const TexWidth = GridDim[0] * Columns;
const TexHeight = Math.ceil(GridDim[2] / Columns) * GridDim[1];
const Resolution = 0.75;
const GridMin = Vec3.create(-1, -2, 0.5);
const GridTransform = Vec4.create(GridMin[0], GridMin[1], GridMin[2], 1 / Resolution);

/** a smooth grid (the layers are), rgba bytes per texel */
function createGrid() {
  const grid = new Uint8Array(TexWidth * TexHeight * 4);
  for (let y = 0; y < TexHeight; ++y) {
    for (let x = 0; x < TexWidth; ++x) {
      const o = (y * TexWidth + x) * 4;
      grid[o] = (x * 11 + y * 7) % 256;
      grid[o + 1] = (x * 5 + y * 13 + 40) % 256;
      grid[o + 2] = 255 - ((x * 3 + y * 9) % 256);
      grid[o + 3] = 128 + ((x + y) % 64);
    }
  }
  return grid;
}

/** the texture's linear filter on texel `(x, y)` of the atlas (texel centers at `.5`, clamped to the edge) */
function filter(grid: Uint8Array, x: number, y: number, c: number) {
  const fx = x - 0.5;
  const fy = y - 0.5;
  const x0 = Math.floor(fx);
  const y0 = Math.floor(fy);
  const tx = fx - x0;
  const ty = fy - y0;
  const at = (xi: number, yi: number) => {
    const cx = Math.min(Math.max(xi, 0), TexWidth - 1);
    const cy = Math.min(Math.max(yi, 0), TexHeight - 1);
    return grid[(cy * TexWidth + cx) * 4 + c] / 255;
  };
  const a = at(x0, y0) * (1 - tx) + at(x0 + 1, y0) * tx;
  const b = at(x0, y0 + 1) * (1 - tx) + at(x0 + 1, y0 + 1) * tx;
  return a * (1 - ty) + b * ty;
}

/**
 * channel `c` at model position `p`, as the mesh shaders sample a layer's grid (`texture3dFrom2dLinear`, its values at
 * the cells' corners: half a texel on in x and y), 0 to 255
 */
function sampleReference(grid: Uint8Array, p: Vec3, c: number) {
  const s = GridTransform[3];
  const pos = [
    ((p[0] - GridMin[0]) * s) / GridDim[0],
    ((p[1] - GridMin[1]) * s) / GridDim[1],
    ((p[2] - GridMin[2]) * s) / GridDim[2],
  ];
  const zSlice0 = Math.floor(pos[2] * GridDim[2]);
  const start0 = zSlice0 * GridDim[0];
  const start1 = start0 + GridDim[0];
  const ox = pos[0] * GridDim[0] + 0.5;
  const oy = pos[1] * GridDim[1] + 0.5;
  const c0 = filter(grid, (start0 % TexWidth) + ox, Math.floor(start0 / TexWidth) * GridDim[1] + oy, c);
  const c1 = filter(grid, (start1 % TexWidth) + ox, Math.floor(start1 / TexWidth) * GridDim[1] + oy, c);
  const t = Math.abs(pos[2] * GridDim[2] - zSlice0);
  return (c0 * (1 - t) + c1 * t) * 255;
}

const VertexCount = 37;
const InstanceCount = 3;
const Translations = [Vec3.create(0, 0, 0), Vec3.create(0.6, -0.4, 0.3), Vec3.create(-0.3, 0.5, 0.9)];

/** positions inside the grid (all instances' too), deterministic */
function createPositions() {
  const positions = new Float32Array(VertexCount * 3);
  for (let i = 0; i < VertexCount; ++i) {
    positions[i * 3] = GridMin[0] + 1 + ((i * 0.37) % 3.2);
    positions[i * 3 + 1] = GridMin[1] + 1 + ((i * 0.53) % 1.6);
    positions[i * 3 + 2] = GridMin[2] + 0.5 + ((i * 0.29) % 2.4);
  }
  return positions;
}

function setup(graphics: GraphicsContext & { compute: NonNullable<GraphicsContext['compute']> }) {
  const { compute } = graphics;
  const grid = createGrid();
  const texture = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
  texture.load({ array: grid, width: TexWidth, height: TexHeight, filter: 'linear' });
  const local = createPositions();
  // two vertices of another mesh before this one's, a matrix of another before its
  const positionData = new Float32Array((VertexCount + 2) * 3);
  positionData.set(local, 6);
  const positions = compute.createBuffer(positionData.byteLength, 'test positions');
  positions.write(positionData);
  const matrices = new Float32Array((InstanceCount + 1) * 16);
  Mat4.toArray(Mat4.fromTranslation(Mat4(), Vec3.create(9, 9, 9)), matrices, 0);
  for (let i = 0; i < InstanceCount; ++i)
    Mat4.toArray(Mat4.fromTranslation(Mat4(), Translations[i]), matrices, (i + 1) * 16);
  const transforms = compute.createBuffer(matrices.byteLength, 'test transforms');
  transforms.write(matrices);
  // the second instance not reached
  const reached = compute.createBuffer(InstanceCount * 4, 'test reached');
  reached.write(new Uint32Array([1, 0, 1]));
  const input = {
    vertexCount: VertexCount,
    instanceCount: InstanceCount,
    positions,
    vertexOffset: 2,
    transforms,
    instanceOffset: 1,
    reached,
    grid: texture,
    gridDim: GridDim,
    gridTexDim: Vec2.create(TexWidth, TexHeight),
    gridTransform: GridTransform,
  };
  const position = (element: number) => {
    const instance = Math.floor(element / VertexCount);
    return Vec3.add(Vec3(), Vec3.fromArray(Vec3(), local, (element % VertexCount) * 3), Translations[instance]);
  };
  /** whether the point's trilinear cells are all in the grid: 0 read outside it */
  const inGrid = (element: number) => {
    const p = position(element);
    return [0, 1, 2].every((k) => {
      const cell = (p[k] - GridMin[k]) * GridTransform[3];
      return cell >= 0 && cell <= GridDim[k] - 1;
    });
  };
  const reference = (element: number, c: number) =>
    Math.floor(element / VertexCount) === 1 || !inGrid(element) ? 0 : sampleReference(grid, position(element), c);
  /**
   * before the grid's last cell in each axis, where the CPU path reads the same texels (the filter clamps at the atlas'
   * edge, the CPU path reads past it)
   */
  const inTile = (element: number) => {
    const p = position(element);
    return [0, 1, 2].every((k) => (p[k] - GridMin[k]) * GridTransform[3] < GridDim[k] - 1);
  };
  // the CPU path's trilinear interpolation between the cells' corners (`getTrilinearlyInterpolated`)
  const transformBuffer = new Float32Array(matrices.subarray(16));
  const cpu = getTrilinearlyInterpolated({
    vertexCount: VertexCount,
    instanceCount: InstanceCount,
    transformBuffer,
    positionBuffer: local,
    colorType: 'volumeInstance',
    grid,
    gridDim: GridDim,
    gridTexDim: Vec2.create(TexWidth, TexHeight),
    gridTransform: GridTransform,
    vertexStride: 3,
    colorStride: 4,
    outputStride: 4,
    instances: new Uint8Array([1, 0, 1]),
  }).array;
  return { input, reference, cpu, inTile, inGrid, destroy: () => [positions, transforms, reached].forEach((b) => b.destroy()) };
}

/** the largest difference to the reference: the filter's weights of the hardware's precision, the bytes rounded */
const Tolerance = 2;

describe.skipIf(!WebGPU.module)('layer sampling per vertex and instance (WebGPU compute)', () => {
  it('samples rgba8 and r8 layers as the mesh shaders sample their grids', async () => {
    const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: 16, height: 16 });
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      if (!isComputeLayerSamplingSupported(graphics)) throw new Error('expected layer sampling support');
      const { compute } = graphics;
      const { input, reference, cpu, inTile, inGrid, destroy } = setup(graphics);
      const count = VertexCount * InstanceCount;
      // some points reach past the grid's last cell: read as 0
      let outside = 0;
      for (let e = 0; e < count; ++e) if (Math.floor(e / VertexCount) !== 1 && !inGrid(e)) outside++;
      expect(outside).toBeGreaterThan(0);

      // rgba8 from element 3: an element a word
      const rgba = compute.createBuffer((count + 3) * 4, 'test rgba8');
      computeLayerSampling(compute, { ...input, itemSize: 4 }, rgba, 3);
      const rgbaBytes = new Uint8Array(await rgba.readAsync());
      let rgbaMax = 0;
      for (let e = 0; e < count; ++e) {
        for (let c = 0; c < 4; ++c) {
          rgbaMax = Math.max(rgbaMax, Math.abs(rgbaBytes[(e + 3) * 4 + c] - reference(e, c)));
        }
      }
      expect(rgbaMax).toBeLessThanOrEqual(Tolerance);
      // as the CPU path samples the grid (the filter between the same corners)
      let cpuMax = 0;
      let compared = 0;
      for (let e = 0; e < count; ++e) {
        if (!inTile(e)) continue;
        compared++;
        for (let c = 0; c < 4; ++c) {
          cpuMax = Math.max(cpuMax, Math.abs(rgbaBytes[(e + 3) * 4 + c] - cpu[e * 4 + c]));
        }
      }
      expect(compared).toBeGreaterThan(count / 2);
      expect(cpuMax).toBeLessThanOrEqual(Tolerance);

      // r8 from element 2: four elements a word, the bytes before and after the range kept
      const words = Math.ceil((count + 2 + 4) / 4);
      const r8 = compute.createBuffer(words * 4, 'test r8');
      r8.write(new Uint8Array(words * 4).fill(0xab));
      computeLayerSampling(compute, { ...input, itemSize: 1 }, r8, 2);
      const r8Bytes = new Uint8Array(await r8.readAsync());
      expect(Array.from(r8Bytes.subarray(0, 2))).toEqual([0xab, 0xab]);
      expect(Array.from(r8Bytes.subarray(count + 2, words * 4)).every((b) => b === 0xab)).toBe(true);
      let r8Max = 0;
      for (let e = 0; e < count; ++e) r8Max = Math.max(r8Max, Math.abs(r8Bytes[e + 2] - reference(e, 0)));
      expect(r8Max).toBeLessThanOrEqual(Tolerance);

      // the unreached instance zero in both
      for (let e = VertexCount; e < 2 * VertexCount; ++e) {
        expect(rgbaBytes[(e + 3) * 4]).toBe(0);
        expect(r8Bytes[e + 2]).toBe(0);
      }
      rgba.destroy();
      r8.destroy();
      destroy();
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('samples data per vertex alone (colors by group: one identity transform) as the CPU path does', async () => {
    const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: 16, height: 16 });
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      if (!isComputeLayerSamplingSupported(graphics)) throw new Error('expected layer sampling support');
      const { compute } = graphics;
      const grid = createGrid();
      const texture = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
      texture.load({ array: grid, width: TexWidth, height: TexHeight, filter: 'linear' });
      const local = createPositions();
      const positions = compute.createBuffer(local.byteLength, 'test positions');
      positions.write(local);
      const transforms = compute.createBuffer(64, 'test transform');
      transforms.write(Mat4.toArray(Mat4.identity(), new Float32Array(16), 0) as Float32Array);
      const reached = compute.createBuffer(4, 'test reached');
      reached.write(new Uint32Array([1]));
      const out = compute.createBuffer(VertexCount * 4, 'test rgba8');
      computeLayerSampling(
        compute,
        {
          vertexCount: VertexCount,
          instanceCount: 1,
          positions,
          vertexOffset: 0,
          transforms,
          instanceOffset: 0,
          reached,
          grid: texture,
          gridDim: GridDim,
          gridTexDim: Vec2.create(TexWidth, TexHeight),
          gridTransform: GridTransform,
          itemSize: 4,
        },
        out,
        0,
      );
      const bytes = new Uint8Array(await out.readAsync());
      // the CPU path's data per vertex (`volume`: no instance, the mesh's own space)
      const cpu = getTrilinearlyInterpolated({
        vertexCount: VertexCount,
        instanceCount: 1,
        transformBuffer: Mat4.toArray(Mat4.identity(), new Float32Array(16), 0) as Float32Array,
        positionBuffer: local,
        colorType: 'volume',
        grid,
        gridDim: GridDim,
        gridTexDim: Vec2.create(TexWidth, TexHeight),
        gridTransform: GridTransform,
        vertexStride: 3,
        colorStride: 4,
        outputStride: 4,
      }).array;
      let refMax = 0;
      let cpuMax = 0;
      for (let v = 0; v < VertexCount; ++v) {
        const p = Vec3.fromArray(Vec3(), local, v * 3);
        for (let c = 0; c < 4; ++c) {
          refMax = Math.max(refMax, Math.abs(bytes[v * 4 + c] - sampleReference(grid, p, c)));
          cpuMax = Math.max(cpuMax, Math.abs(bytes[v * 4 + c] - cpu[v * 4 + c]));
        }
      }
      expect(refMax).toBeLessThanOrEqual(Tolerance);
      expect(cpuMax).toBeLessThanOrEqual(Tolerance);
      expect(bytes.some((b) => b > 0)).toBe(true);
      for (const b of [positions, transforms, reached, out]) b.destroy();
      texture.destroy();
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
