/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import { estimateIsosurfaceVertexCount } from '@molstar/graphics/geo/util/marching-cubes/estimate';
import type { ComputeContext } from '@molstar/graphics/gl/compute';
import {
  IsosurfaceChunkSize,
  IsosurfaceTableHeader,
  IsosurfaceVertexBytes,
} from '@molstar/graphics/gl/compute/marching-cubes/chunk-size';
import {
  type ChunkedSurface,
  computeChunkedIsosurface,
  destroyChunkedSurface,
  getChunkPages,
  getVertexArenas,
  growVertexArena,
  readChunkedVertices,
  reserveVertexArenaPages,
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
import { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';
import { isNullTexture } from '@molstar/graphics/gl/texture';
import { createWrappedTensor } from '@molstar/graphics/repr/volume/util';
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

/** A periodic density, `n` minus the distance to `c` across the faces: a sphere of radius `r` at iso value `n - r`, cut by
 * the faces near `c` and continued on the opposite ones */
function periodicSphereValues(n: number, c: number) {
  const values = new Float32Array(n * n * n);
  const d = (v: number) => Math.min(Math.abs(v - c), n - Math.abs(v - c));
  for (let z = 0; z < n; ++z) {
    for (let y = 0; y < n; ++y) {
      for (let x = 0; x < n; ++x) values[x + y * n + z * n * n] = n - Math.hypot(d(x), d(y), d(z));
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

/** The positions, normals and groups of a surface's first `count` vertices, read through its page table */
async function readVertices(compute: ComputeContext, surface: ChunkedSurface, count: number) {
  const chunks = { table: surface.table, vertices: surface.arena.vertices.ref.value };
  const { vertices, normals, groups, vertexCount } = await readChunkedVertices(compute, chunks, count);
  const out: { p: Vec3; n: Vec3; group: number }[] = [];
  for (let i = 0; i < vertexCount; ++i) {
    out.push({
      p: Vec3.fromArray(Vec3(), vertices, i * 4),
      n: Vec3.fromArray(Vec3(), normals, i * 4),
      // as `unpackRGBToInt`
      group: (groups[i * 4] << 16) + (groups[i * 4 + 1] << 8) + groups[i * 4 + 2] - 1,
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
      // indexed, a vertex per crossed edge as the CPU marching cubes
      expect(state.vertices).toBe(cpu.vertexCount);
      expect(state.pages).toBe(getChunkPages(state.vertices) + getChunkPages(state.count / 3));
      expect(state.missingPages).toBe(0);
      const c = (n - 1) / 2;
      for (const { p, n: normal, group } of await readVertices(compute, surface, state.count)) {
        // a cell of the grid (the cell index of buffer inputs)
        expect(group).toBeGreaterThanOrEqual(0);
        expect(group).toBeLessThan(n * n * n);
        const radial = Vec3.sub(Vec3(), p, Vec3.create(c, c, c));
        expect(Math.abs(Vec3.magnitude(radial) - r)).toBeLessThan(0.1);
        expect(Math.abs(Vec3.magnitude(normal) - 1)).toBeLessThan(0.01);
        // out of the sphere, against the density gradient (as the CPU normals): radial, within the 24-bit octahedral
        // encoding and the gradient's error
        expect(Vec3.dot(normal, Vec3.normalize(radial, radial))).toBeGreaterThan(0.995);
      }
    });
  });

  it('wraps the grid of a periodic volume around its data, as the CPU marching cubes', async () => {
    await withCompute(async (graphics) => {
      const { compute } = graphics;
      const n = 24;
      const c = 1.5;
      const r = 6.25;
      const values = periodicSphereValues(n, c);
      const volume = compute.createBuffer(values.byteLength, 'test volume');
      volume.write(values);
      const dataDim = Vec3.create(n, n, n);
      const input: IsosurfaceInput = { kind: 'buffer', volume, strides: Vec3.create(1, n, n * n), dataDim };
      const tensor = Tensor.create(Tensor.Space([n, n, n], [2, 1, 0], Float32Array), Tensor.Data1(values));
      const cpu = await computeMarchingCubesMesh({ isoLevel: n - r, scalarField: createWrappedTensor(tensor) }).run();
      // a point more per axis: the cells across the faces
      const gridDim = Vec3.create(n + 1, n + 1, n + 1);
      const { surface, count } = await computeChunkedIsosurface(compute, input, gridDim, n - r, Mat4.identity(), false);
      const state = await count;
      expect(state.count).toBe(cpu.triangleCount * 3);
      // on the sphere around `c` or one of its periodic images
      const d = (v: number) => Math.min(Math.abs(v - c), n - Math.abs(v - c));
      for (const { p } of await readVertices(compute, surface, state.count)) {
        expect(Math.abs(Math.hypot(d(p[0]), d(p[1]), d(p[2])) - r)).toBeLessThan(0.1);
      }
      volume.destroy();
    });
  });

  it('fills a texture mesh last filled into textures (another context) into the arena, without its textures', async () => {
    await withCompute(async (graphics) => {
      const { compute } = graphics;
      const n = 16;
      const texture = graphics.resources.texture('image-float32', 'rgba', 'float', 'nearest');
      texture.define(4, 4);
      const mesh = TextureMesh.create(3, 1, texture, texture, texture, Sphere3D(), undefined);
      const { surface, count } = await computeChunkedIsosurface(
        compute,
        sphereVolume(compute, n),
        Vec3.create(n, n, n),
        n - 4,
        Mat4.identity(),
        false,
      );
      const state = await count;
      TextureMesh.createChunked(state.count, 1, surface, compute, Sphere3D(), mesh);
      expect(isNullTexture(mesh.vertexTexture.ref.value)).toBe(true);
      expect(isNullTexture(mesh.groupTexture.ref.value)).toBe(true);
      expect(isNullTexture(mesh.normalTexture.ref.value)).toBe(true);
      expect(mesh.isoChunks.ref.value).toBe(true);
      expect(mesh.doubleBuffer.get(graphics)).toBeUndefined();
      mesh.doubleBuffer.destroy();
      texture.destroy();
    });
  });

  it('reuses its pages on updates and returns them when destroyed', async () => {
    await withCompute(async (graphics) => {
      const { compute } = graphics;
      const n = 40;
      const input = sphereVolume(compute, n);
      const dim = Vec3.create(n, n, n);
      const first = await computeChunkedIsosurface(compute, input, dim, n - 15, Mat4.identity(), false);
      const { arena } = first.surface;
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
      // without estimates or reservations, the arena isn't grown beforehand: of two spheres of about 36 pages (22 K vertices, 44 K triangles), the
      // second doesn't fit the 64 pages of the first arena
      const iso = n - 50;
      const a = await computeChunkedIsosurface(compute, input, dim, iso, Mat4.identity(), false);
      const b = await computeChunkedIsosurface(compute, input, dim, iso, Mat4.identity(), false);
      const [stateA, stateB] = await Promise.all([a.count, b.count]);
      expect(stateA.missingPages).toBe(0);
      expect(stateB.missingPages).toBeGreaterThan(0);
      expect(stateB.drawn).toBeLessThan(stateB.count);
      expect(stateB.drawn % 3).toBe(0);

      expect(b.surface.arena).toBe(getVertexArenas(compute)[0]);
      expect(await growVertexArena(compute, b.surface.arena, stateB.missingPages)).toBe(true);
      const again = await computeChunkedIsosurface(compute, input, dim, iso, Mat4.identity(), false, b.surface);
      const stateAgain = await again.count;
      expect(stateAgain.missingPages).toBe(0);
      expect(stateAgain.drawn).toBe(stateA.count);
      // the first surface's vertices survived the growth (copied into the new buffers)
      const vertices = await readVertices(compute, a.surface, stateA.count);
      const c = (n - 1) / 2;
      for (const { p } of vertices) expect(Math.abs(Vec3.distance(p, Vec3.create(c, c, c)) - 50)).toBeLessThan(0.1);
    });
  });

  it('gets all pages at once with the pages of an estimate from the grid reserved', async () => {
    await withCompute(async (graphics) => {
      const { compute } = graphics;
      const n = 128;
      const input = sphereVolume(compute, n);
      const dim = Vec3.create(n, n, n);
      const iso = n - 50;
      // the two spheres of the case above, the arena grown beforehand by the estimate (sampled: at least the count)
      const estimate = estimateIsosurfaceVertexCount(sphereValues(n), dim, [1, n, n * n], iso);
      const extract = async () => {
        const reservation = await reserveVertexArenaPages(compute, getChunkPages(estimate));
        return computeChunkedIsosurface(
          compute,
          input,
          dim,
          iso,
          Mat4.identity(),
          false,
          undefined,
          estimate,
          undefined,
          reservation,
        );
      };
      const a = await extract();
      const b = await extract();
      const [stateA, stateB] = await Promise.all([a.count, b.count]);
      expect(estimate).toBeGreaterThanOrEqual(stateA.count);
      expect(stateA.missingPages).toBe(0);
      expect(stateB.missingPages).toBe(0);
      expect(stateB.drawn).toBe(stateB.count);
    });
  });

  it('extracts a surface its estimate fell short of again, with a table sized to its count', async () => {
    await withCompute(async (graphics) => {
      const { compute } = graphics;
      const n = 128;
      const input = sphereVolume(compute, n);
      const dim = Vec3.create(n, n, n);
      const extract = (target: Parameters<typeof extractIsosurface>[6]) =>
        extractIsosurface(graphics, input, dim, n - 50, Mat4.identity(), false, target);
      const sphere = Sphere3D.create(Vec3.create(n / 2, n / 2, n / 2), n);
      // a table of 18 pages for the estimate of one vertex, the sphere about 45
      const mesh = await computeTextureMeshIsosurface(graphics, extract, 1, sphere, undefined, 1);
      const firstTable = mesh.isoTable.ref.value;
      await mesh.meta.extracted;
      const { table } = mesh.doubleBuffer.getChunks(graphics)!;
      // 18 vertex and 18 triangle pages
      expect(firstTable.byteLength).toBe((IsosurfaceTableHeader + 18 + 18) * 4);
      expect(mesh.isoTable.ref.value).toBe(table);
      const header = new Uint32Array(await table.readAsync(0, 16));
      // no overflow, all drawn
      expect(header[2]).toBe(0);
      expect(mesh.meta.readVertexCount).toBe(header[3]);
      mesh.doubleBuffer.destroy();
    });
  });

  /** `graphics` with arenas of at most the 64 pages of the first (its growth beyond them a new arena), `oom`: no memory */
  function limitArenas(graphics: GraphicsContext & { compute: ComputeContext }, oom: boolean) {
    const compute = Object.create(graphics.compute, {
      maxStorageBufferSize: { value: 64 * IsosurfaceChunkSize * IsosurfaceVertexBytes },
      ...(oom ? { tryCreateBuffer: { value: async () => undefined } } : {}),
    }) as ComputeContext;
    return Object.create(graphics, { compute: { value: compute } }) as typeof graphics;
  }

  it('moves a surface its arena lacks the pages for into a new arena, without copying', async () => {
    await withCompute(async (graphics) => {
      const limited = limitArenas(graphics, false);
      const n = 128;
      const input = sphereVolume(limited.compute, n);
      const dim = Vec3.create(n, n, n);
      const extract = (target: Parameters<typeof extractIsosurface>[6]) =>
        extractIsosurface(limited, input, dim, n - 50, Mat4.identity(), false, target);
      const sphere = Sphere3D.create(Vec3.create(n / 2, n / 2, n / 2), n);
      // two spheres of about 36 pages (22 K vertices, 44 K triangles): the second doesn't fit the rest of the first arena's 64, extracted again in
      // another arena
      const a = await computeTextureMeshIsosurface(limited, extract, 1, sphere);
      const b = await computeTextureMeshIsosurface(limited, extract, 1, sphere);
      await Promise.all([a.meta.extracted, b.meta.extracted]);
      const arenas = getVertexArenas(limited.compute);
      expect(arenas.length).toBe(2);
      const surfaceA = a.doubleBuffer.getChunks(limited)!;
      const surfaceB = b.doubleBuffer.getChunks(limited)!;
      expect(surfaceA.arena).toBe(arenas[0]);
      expect(surfaceB.arena).toBe(arenas[1]);
      // what draws it binds the new arena's vertices
      expect(b.isoVertices.ref.value).toBe(arenas[1].vertices.ref.value);
      expect(b.meta.readVertexCount).toBe(a.meta.readVertexCount);
      // its pages returned to the first arena
      expect(new Uint32Array(await arenas[0].free.readAsync(0, 4))[0]).toBe(64 - surfaceA.pages);
      const c = (n - 1) / 2;
      for (const { p } of await readVertices(limited.compute, surfaceB, b.meta.readVertexCount!)) {
        expect(Math.abs(Vec3.distance(p, Vec3.create(c, c, c)) - 50)).toBeLessThan(0.1);
      }
      a.doubleBuffer.destroy();
      b.doubleBuffer.destroy();
    });
  });

  it('keeps its arenas while a surface is created again every frame (a trajectory played back)', async () => {
    await withCompute(async (graphics) => {
      const limited = limitArenas(graphics, false);
      const n = 128;
      const input = sphereVolume(limited.compute, n);
      const dim = Vec3.create(n, n, n);
      const extract = (target: Parameters<typeof extractIsosurface>[6]) =>
        extractIsosurface(limited, input, dim, n - 50, Mat4.identity(), false, target);
      const sphere = Sphere3D.create(Vec3.create(n / 2, n / 2, n / 2), n);
      // the new frame's surface before the old one is destroyed: two at once, more than an arena of 64 pages holds
      let previous = await computeTextureMeshIsosurface(limited, extract, 1, sphere);
      await previous.meta.extracted;
      let arenas: readonly unknown[] | undefined;
      for (let frame = 0; frame < 5; ++frame) {
        const next = await computeTextureMeshIsosurface(limited, extract, 1, sphere);
        await next.meta.extracted;
        previous.doubleBuffer.destroy();
        previous = next;
        expect(next.meta.readVertexCount).toBeGreaterThan(0);
        // the same arenas from the first frame on: none released and created again
        const current = [...getVertexArenas(limited.compute)];
        if (arenas) expect(current).toEqual(arenas);
        arenas = current;
      }
      expect(arenas).toHaveLength(2);
      previous.doubleBuffer.destroy();
    });
  });

  it('reports surfaces drawn incomplete when the device has no memory for another arena', async () => {
    await withCompute(async (graphics) => {
      const limited = limitArenas(graphics, true);
      const n = 128;
      const input = sphereVolume(limited.compute, n);
      const dim = Vec3.create(n, n, n);
      const extract = (target: Parameters<typeof extractIsosurface>[6]) =>
        extractIsosurface(limited, input, dim, n - 50, Mat4.identity(), false, target);
      const sphere = Sphere3D.create(Vec3.create(n / 2, n / 2, n / 2), n);
      const reported = new Promise<string>((resolve) => {
        const sub = graphics.resources.outOfMemory.subscribe((m) => {
          sub.unsubscribe();
          resolve(m);
        });
      });
      // two spheres of about 36 pages (22 K vertices, 44 K triangles): the second gets the rest of the 64
      const a = await computeTextureMeshIsosurface(limited, extract, 1, sphere);
      const b = await computeTextureMeshIsosurface(limited, extract, 1, sphere);
      await Promise.all([a.meta.extracted, b.meta.extracted]);
      expect(getVertexArenas(limited.compute).length).toBe(1);
      expect(b.meta.readVertexCount).toBeLessThan(a.meta.readVertexCount!);
      expect(await reported).toMatch(/^Not enough GPU memory: 1 surface drawn incomplete.*the GPU has no memory/);
      a.doubleBuffer.destroy();
      b.doubleBuffer.destroy();
    });
  });

  it('binds the grown page table of a surface extracted again with its count', async () => {
    await withCompute(async (graphics) => {
      const { compute } = graphics;
      const n = 128;
      const input = sphereVolume(compute, n);
      const dim = Vec3.create(n, n, n);
      const extract = (target: Parameters<typeof extractIsosurface>[6]) =>
        extractIsosurface(graphics, input, dim, n - 50, Mat4.identity(), false, target);
      const sphere = Sphere3D.create(Vec3.create(n / 2, n / 2, n / 2), n);
      // without estimates the second mesh doesn't get all its pages, with an estimate of one vertex its table is too
      // small: extracted again with its count, after the arena grew, into a larger table
      const a = await computeTextureMeshIsosurface(graphics, extract, 1, sphere);
      const b = await computeTextureMeshIsosurface(graphics, extract, 1, sphere, undefined, 1);
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
