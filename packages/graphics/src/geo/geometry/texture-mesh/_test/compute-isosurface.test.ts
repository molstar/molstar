/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Subject } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import type { ComputeContext, StorageBuffer } from '@molstar/graphics/gl/compute';
import { IsosurfaceChunkSize, IsosurfaceVertexBytes } from '@molstar/graphics/gl/compute/marching-cubes/chunk-size';
import {
  type ChunkedSurface,
  type ChunkedSurfaceCount,
  getVertexArenas,
} from '@molstar/graphics/gl/compute/marching-cubes/chunks';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { computeTextureMeshIsosurface, type IsosurfaceExtraction } from '../compute-isosurface.js';

function fakeBuffer(byteLength = 16): StorageBuffer {
  return {
    id: 0,
    byteLength,
    write: () => {},
    clear: () => {},
    copyTo: () => {},
    readAsync: async () => new ArrayBuffer(byteLength),
    getByteCount: () => byteLength,
    destroy: () => {},
  };
}

/**
 * A compute context for the vertex arena's bookkeeping: buffers and kernels that do nothing; `maxStorageBufferSize`
 * bounds the arena's pages (64, of a new arena, at 3 MiB)
 */
function fakeGraphics(maxStorageBufferSize = 1 << 30) {
  const compute = {
    maxStorageBufferSize,
    createBuffer: (byteLength: number) => fakeBuffer(byteLength),
    tryCreateBuffer: async (byteLength: number) => fakeBuffer(byteLength),
    createKernel: () => ({ dispatch: () => {}, destroy: () => {} }),
    beginComputePass: () => {},
    endComputePass: () => {},
  } as unknown as ComputeContext;
  return {
    compute,
    resources: { pending: { track: () => {} }, outOfMemory: new Subject<string>() },
  } as unknown as GraphicsContext & { compute: ComputeContext };
}

/**
 * Extractions as `computeChunkedIsosurface` does them: into the surface given (its table and draw arguments reused, as
 * for a mesh extracted again), else a new one; each extraction's count is read back when `resolve` is called for it
 */
function fakeExtractions(compute: ComputeContext) {
  const pending: ((count: ChunkedSurfaceCount) => void)[] = [];
  const extract: IsosurfaceExtraction = async ({ surface }) => {
    const [arena] = getVertexArenas(compute);
    const s = (surface ?? {
      arena,
      table: fakeBuffer(),
      retired: [],
      drawArgs: fakeBuffer(),
      pages: 0,
    }) as ChunkedSurface;
    return {
      surface: s,
      drawIndirect: { buffer: s.drawArgs, byteOffset: 0 },
      count: new Promise<ChunkedSurfaceCount>((r) => pending.push(r)),
    };
  };
  /** the count of the `i`th extraction, `missingPages` of its pages not available (`arenaFull`: for want of free pages) */
  const resolve = async (i: number, drawn: number, missingPages = 0, arenaFull = false) => {
    const vertices = Math.ceil(drawn / 6);
    pending[i]({ count: drawn, vertices, drawn: missingPages ? 0 : drawn, pages: 1, missingPages, arenaFull });
    // the readback's handling: reservations in turn (a promise chain) before an extraction again
    for (let k = 0; k < 20; ++k) await Promise.resolve();
  };
  return { extract, resolve, calls: () => pending.length };
}

// GPU surfaces are drawn with the count on the GPU until it is read back, then directly with it: indirect draws cost
// much more GPU time in Chrome, over hundreds of surfaces most of a frame (`TextureMesh.setDrawnVertexCount`)
describe('compute isosurface draws', () => {
  it('draws indirectly until the vertex count is read back, then directly with it', async () => {
    const graphics = fakeGraphics();
    const x = fakeExtractions(graphics.compute);
    const mesh = await computeTextureMeshIsosurface(graphics, x.extract, 1, Sphere3D(), undefined, 9000);
    const args = mesh.drawIndirect.ref.value!;
    expect(args.count).toBeUndefined();
    expect(mesh.vertexCount).toBe(9000);

    await x.resolve(0, 4242);
    expect(mesh.drawIndirect.ref.value?.buffer).toBe(args.buffer);
    expect(mesh.drawIndirect.ref.value?.count).toBe(4242);
    // the upper bound stays: the stride of vertex-instance data laid out with it
    expect(mesh.vertexCount).toBe(9000);
    expect(mesh.vertexCountCell.ref.value).toBe(9000);
  });

  it('draws indirectly again when the mesh is extracted again into the same draw arguments', async () => {
    const graphics = fakeGraphics();
    const x = fakeExtractions(graphics.compute);
    const mesh = await computeTextureMeshIsosurface(graphics, x.extract, 1, Sphere3D(), undefined, 9000);
    const buffer = mesh.drawIndirect.ref.value!.buffer;
    await x.resolve(0, 4242);

    const again = await computeTextureMeshIsosurface(graphics, x.extract, 1, Sphere3D(), mesh);
    expect(again).toBe(mesh);
    // the same surface and arguments buffer, the count of the previous extraction dropped
    expect(mesh.drawIndirect.ref.value?.buffer).toBe(buffer);
    expect(mesh.drawIndirect.ref.value?.count).toBeUndefined();
    await x.resolve(1, 5001);
    expect(mesh.drawIndirect.ref.value?.count).toBe(5001);
  });

  it('ignores the count of an extraction a newer one replaced', async () => {
    const graphics = fakeGraphics();
    const x = fakeExtractions(graphics.compute);
    const mesh = await computeTextureMeshIsosurface(graphics, x.extract, 1, Sphere3D(), undefined, 9000);
    await computeTextureMeshIsosurface(graphics, x.extract, 1, Sphere3D(), mesh, 9000);
    // the first count arrives after the second extraction started into the same surface: only the extraction is newer
    await x.resolve(0, 100);
    expect(mesh.drawIndirect.ref.value?.count).toBeUndefined();
    await x.resolve(1, 300);
    expect(mesh.drawIndirect.ref.value?.count).toBe(300);
  });

  it('draws the count of the extraction run again after the arena grew, not of the incomplete one', async () => {
    const graphics = fakeGraphics();
    const x = fakeExtractions(graphics.compute);
    const mesh = await computeTextureMeshIsosurface(graphics, x.extract, 1, Sphere3D(), undefined, 3000);
    await x.resolve(0, 12_000, 2);
    expect(x.calls()).toBe(2);
    expect(mesh.drawIndirect.ref.value?.count).toBeUndefined();
    await x.resolve(1, 12_000);
    expect(mesh.drawIndirect.ref.value?.count).toBe(12_000);
    // beyond the estimate: the upper bound grew to it
    expect(mesh.vertexCount).toBe(12_000);
  });

  it('draws a surface the arena at its limit lacks pages for as far as it got, and reports it once', async () => {
    vi.useFakeTimers();
    // the 64 pages of a new arena, the most the storage buffer allows
    const graphics = fakeGraphics(64 * IsosurfaceChunkSize * IsosurfaceVertexBytes);
    const reports: string[] = [];
    graphics.resources.outOfMemory.subscribe((m) => reports.push(m));
    const x = fakeExtractions(graphics.compute);
    const a = await computeTextureMeshIsosurface(graphics, x.extract, 1, Sphere3D(), undefined, 3000);
    const b = await computeTextureMeshIsosurface(graphics, x.extract, 1, Sphere3D(), undefined, 3000);
    await x.resolve(0, 400_000, 100, true);
    await x.resolve(1, 400_000, 100, true);
    // not extracted again: the arena can't grow and has fewer free pages than missing
    expect(x.calls()).toBe(2);
    for (const mesh of [a, b]) expect(mesh.drawIndirect.ref.value?.count).toBe(0);
    expect(reports).toEqual([]);
    await vi.runAllTimersAsync();
    expect(reports).toHaveLength(1);
    expect(reports[0]).toMatch(/^Not enough GPU memory: 2 surfaces drawn incomplete/);
  });

  it('draws an empty surface directly, with no vertices', async () => {
    const graphics = fakeGraphics();
    const x = fakeExtractions(graphics.compute);
    const mesh = await computeTextureMeshIsosurface(graphics, x.extract, 1, Sphere3D(), undefined, 3000);
    await x.resolve(0, 0);
    expect(mesh.drawIndirect.ref.value?.count).toBe(0);
  });
});
