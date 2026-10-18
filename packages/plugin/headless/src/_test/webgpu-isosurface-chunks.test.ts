/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import type { ComputeContext } from '@molstar/graphics/gl/compute';
import { IsosurfaceChunkSize } from '@molstar/graphics/gl/compute/marching-cubes/chunk-size';
import {
  type ChunkedSurface,
  computeChunkedIsosurface,
  destroyChunkedSurface,
  getVertexArena,
  growVertexArena,
} from '@molstar/graphics/gl/compute/marching-cubes/chunks';
import type { IsosurfaceInput } from '@molstar/graphics/gl/compute/marching-cubes/compute';
import {
  computeTextureMeshIsosurface,
  extractIsosurface,
} from '@molstar/graphics/geo/geometry/texture-mesh/compute-isosurface';
import { Sphere3D } from '@molstar/core/math/geometry';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

/** A density falling off from the grid center, `n` minus the distance: a sphere of radius `r` at iso value `n - r` */
function sphereValues(n: number) {
  const values = new Float32Array(n * n * n);
  const c = (n - 1) / 2;
  for (let z = 0; z < n; ++z) {
    for (let y = 0; y < n; ++y) {
      for (let x = 0; x < n; ++x) values[x + y * n + z * n * n] = n - Math.hypot(x - c, y - c, z - c);
    }
  }
  return values;
}

function sphereVolume(compute: ComputeContext, n: number): IsosurfaceInput {
  const values = sphereValues(n);
  const volume = compute.createBuffer(values.byteLength, 'test volume');
  volume.write(values);
  return { kind: 'buffer', volume, strides: Vec3.create(1, n, n * n) };
}

/** The positions and normals of a surface's first `count` vertices, read through its page table */
async function readVertices(compute: ComputeContext, surface: ChunkedSurface, count: number) {
  const arena = getVertexArena(compute);
  const table = new Uint32Array(await surface.table.readAsync());
  const positions = new Float32Array(await arena.positions.ref.value.readAsync());
  const normals = new Uint32Array(await arena.normals.ref.value.readAsync());
  const snorm = (v: number) => Math.max(-1, ((v << 16) >> 16) / 32767);
  const out: { p: Vec3; n: Vec3 }[] = [];
  for (let i = 0; i < count; ++i) {
    const slot = table[4 + Math.floor(i / IsosurfaceChunkSize)] * IsosurfaceChunkSize + (i % IsosurfaceChunkSize);
    const [nxy, nz] = [normals[slot * 2], normals[slot * 2 + 1]];
    out.push({
      p: Vec3.create(positions[slot * 4], positions[slot * 4 + 1], positions[slot * 4 + 2]),
      n: Vec3.create(snorm(nxy & 0xffff), snorm(nxy >>> 16), snorm(nz & 0xffff)),
    });
  }
  return out;
}

async function withCompute(fn: (graphics: GraphicsContext & { compute: ComputeContext }) => Promise<void>) {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: 16, height: 16 },
    { backend: 'webgpu' },
  );
  try {
    await plugin.init();
    const graphics = plugin.canvas3d!.graphics;
    if (!graphics.compute) throw new Error('no compute');
    await fn(graphics as GraphicsContext & { compute: ComputeContext });
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  } finally {
    plugin.dispose();
  }
}

describe.skipIf(!WebGPU.module)('chunked isosurfaces on WebGPU', () => {
  it('extracts into pages of the arena, as many vertices as the CPU marching cubes, on the surface', async () => {
    await withCompute(async (graphics) => {
      const { compute } = graphics;
      const n = 32;
      const r = 10;
      const input = sphereVolume(compute, n);
      const iso = n - r;
      // x fastest, as the strides of the storage buffer
      const scalarField = Tensor.create(
        Tensor.Space([n, n, n], [2, 1, 0], Float32Array),
        Tensor.Data1(sphereValues(n)),
      );
      const cpu = await computeMarchingCubesMesh({ isoLevel: iso, scalarField }).run();
      const { surface, count } = await computeChunkedIsosurface(
        compute,
        input,
        Vec3.create(n, n, n),
        iso,
        Mat4.identity(),
        false,
      );
      const state = await count;
      expect(state.count).toBe(cpu.triangleCount * 3);
      expect(state.drawn).toBe(state.count);
      expect(state.pages).toBe(Math.ceil(state.count / IsosurfaceChunkSize));
      expect(state.missingPages).toBe(0);
      const c = (n - 1) / 2;
      for (const { p, n: normal } of await readVertices(compute, surface, state.count)) {
        const radial = Vec3.sub(Vec3(), p, Vec3.create(c, c, c));
        expect(Math.abs(Vec3.magnitude(radial) - r)).toBeLessThan(0.1);
        expect(Math.abs(Vec3.magnitude(normal) - 1)).toBeLessThan(0.01);
        // out of the sphere, against the density gradient (as the CPU normals)
        expect(Vec3.dot(normal, radial)).toBeGreaterThan(0);
      }
    });
  });

  it('reuses its pages on updates and returns them when destroyed', async () => {
    await withCompute(async (graphics) => {
      const { compute } = graphics;
      const n = 40;
      const input = sphereVolume(compute, n);
      const dim = Vec3.create(n, n, n);
      const arena = getVertexArena(compute);
      const first = await computeChunkedIsosurface(compute, input, dim, n - 15, Mat4.identity(), false);
      const large = await first.count;
      const second = await computeChunkedIsosurface(compute, input, dim, n - 6, Mat4.identity(), false, first.surface);
      const small = await second.count;
      expect(second.surface).toBe(first.surface);
      expect(small.pages).toBeLessThan(large.pages);
      expect(arena.freeRead).toBe(arena.pageCount - small.pages);
      destroyChunkedSurface(compute, second.surface);
      const free = new Uint32Array(await arena.free.readAsync(0, 4))[0];
      expect(free).toBe(arena.pageCount);
    });
  });

  it('reports pages a surface did not get, and gets them after the arena grew', async () => {
    await withCompute(async (graphics) => {
      const { compute } = graphics;
      const n = 128;
      const input = sphereVolume(compute, n);
      const dim = Vec3.create(n, n, n);
      const arena = getVertexArena(compute);
      // without estimates, the arena isn't grown beforehand: of two spheres of about 45 pages, the second doesn't fit
      // the 64 pages of a new arena
      const iso = n - 40;
      const a = await computeChunkedIsosurface(compute, input, dim, iso, Mat4.identity(), false);
      const b = await computeChunkedIsosurface(compute, input, dim, iso, Mat4.identity(), false);
      const [stateA, stateB] = await Promise.all([a.count, b.count]);
      expect(stateA.missingPages).toBe(0);
      expect(stateB.missingPages).toBeGreaterThan(0);
      expect(stateB.drawn).toBeLessThan(stateB.count);
      expect(stateB.drawn % 3).toBe(0);

      expect(growVertexArena(compute, arena, stateB.missingPages)).toBe(true);
      const again = await computeChunkedIsosurface(compute, input, dim, iso, Mat4.identity(), false, b.surface);
      const stateAgain = await again.count;
      expect(stateAgain.missingPages).toBe(0);
      expect(stateAgain.drawn).toBe(stateA.count);
      // the first surface's vertices survived the growth (copied into the new buffers)
      const vertices = await readVertices(compute, a.surface, stateA.count);
      const c = (n - 1) / 2;
      for (const { p } of vertices) expect(Math.abs(Vec3.distance(p, Vec3.create(c, c, c)) - 40)).toBeLessThan(0.1);
    });
  });

  it('binds the grown page table of a surface extracted again after the arena grew', async () => {
    await withCompute(async (graphics) => {
      const { compute } = graphics;
      const n = 128;
      const input = sphereVolume(compute, n);
      const dim = Vec3.create(n, n, n);
      const extract = (target: Parameters<typeof extractIsosurface>[6]) =>
        extractIsosurface(graphics, input, dim, n - 40, Mat4.identity(), false, target);
      const sphere = Sphere3D.create(Vec3.create(n / 2, n / 2, n / 2), n);
      // without estimates the second mesh doesn't get all its pages: extracted again after the arena grew, its table
      // grows with the arena
      const a = await computeTextureMeshIsosurface(graphics, extract, 1, sphere);
      const b = await computeTextureMeshIsosurface(graphics, extract, 1, sphere);
      const firstTable = b.isoTable.ref.value;
      await Promise.all([a.meta.extracted, b.meta.extracted]);
      expect(b.doubleBuffer.getChunks(graphics)!.table).not.toBe(firstTable);
      // what draws the meshes binds their current tables (the tables grown out of are released)
      for (const mesh of [a, b]) expect(mesh.isoTable.ref.value).toBe(mesh.doubleBuffer.getChunks(graphics)!.table);
      for (const mesh of [a, b]) expect(mesh.doubleBuffer.getChunks(graphics)!.retired).toEqual([]);
      a.doubleBuffer.destroy();
      b.doubleBuffer.destroy();
    });
  });
});
