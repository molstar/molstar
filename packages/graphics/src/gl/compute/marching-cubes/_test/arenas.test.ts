/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { ComputeContext, StorageBuffer } from '../../../compute.js';
import {
  computeChunkedIsosurface,
  destroyChunkedSurface,
  getVertexArenas,
  releaseVertexArenaPages,
  reserveVertexArenaPages,
} from '../chunks.js';

function fakeBuffer(byteLength = 16): StorageBuffer & { destroyed: boolean } {
  const buffer = {
    id: 0,
    byteLength,
    destroyed: false,
    write: () => {},
    clear: () => {},
    copyTo: () => {},
    readAsync: async () => new ArrayBuffer(byteLength),
    getByteCount: () => byteLength,
    destroy: () => {
      buffer.destroyed = true;
    },
  };
  return buffer;
}

/** A compute context for the arenas' bookkeeping: buffers and kernels that do nothing, read back as zeros */
function fakeCompute() {
  return {
    maxStorageBufferSize: 1 << 30,
    createBuffer: (byteLength: number) => fakeBuffer(byteLength),
    tryCreateBuffer: async (byteLength: number) => fakeBuffer(byteLength),
    createKernel: () => ({ dispatch: () => {}, destroy: () => {} }),
    hasKernels: () => true,
    beginComputePass: () => {},
    endComputePass: () => {},
  } as unknown as ComputeContext;
}

/** More pages than an arena of 256 MiB has: an arena of its own */
const LargePages = 6000;

const extract = (compute: ComputeContext, reservation: Awaited<ReturnType<typeof reserveVertexArenaPages>>) => {
  const input = { kind: 'buffer', volume: fakeBuffer(4), strides: Vec3.create(1, 1, 1) } as const;
  const dim = Vec3.create(2, 2, 2);
  return computeChunkedIsosurface(compute, input, dim, 0, Mat4.identity(), false, undefined, 1, undefined, reservation);
};

// Arenas beyond the first are 256 MiB or more: released once unused, but not between the frames of a trajectory
// played back, whose surfaces may be created again each frame
describe('vertex arenas', () => {
  it('releases an arena without reservations or surfaces after a while, keeps it when used again in between', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout'] });
    const compute = fakeCompute();
    const first = await reserveVertexArenaPages(compute, LargePages);
    const [, arena] = getVertexArenas(compute);
    expect(first.arena).toBe(arena);
    releaseVertexArenaPages(first);
    vi.advanceTimersByTime(2900);
    expect(getVertexArenas(compute)).toContain(arena);

    // the next frame's surface takes the pages again
    const again = await reserveVertexArenaPages(compute, LargePages);
    expect(again.arena).toBe(arena);
    vi.advanceTimersByTime(10_000);
    expect(getVertexArenas(compute)).toContain(arena);

    releaseVertexArenaPages(again);
    vi.advanceTimersByTime(3000);
    expect(getVertexArenas(compute)).not.toContain(arena);
    expect((arena.vertices.ref.value as ReturnType<typeof fakeBuffer>).destroyed).toBe(true);
  });

  it('keeps an arena while a surface has pages there, releases it after the surface is destroyed', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout'] });
    const compute = fakeCompute();
    const result = await extract(compute, await reserveVertexArenaPages(compute, LargePages));
    const { arena } = result.surface;
    // the reservation returned with the count read back, the surface's pages still there
    await result.count;
    vi.advanceTimersByTime(10_000);
    expect(getVertexArenas(compute)).toContain(arena);

    destroyChunkedSurface(compute, result.surface);
    vi.advanceTimersByTime(3000);
    expect(getVertexArenas(compute)).not.toContain(arena);
  });

  it('creates the first arena again after all were released', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout'] });
    const compute = fakeCompute();
    const result = await extract(compute, await reserveVertexArenaPages(compute, 1));
    const [first] = getVertexArenas(compute);
    expect(result.surface.arena).toBe(first);
    await result.count;
    destroyChunkedSurface(compute, result.surface);
    vi.advanceTimersByTime(3000);
    const arenas = getVertexArenas(compute);
    expect(arenas).toHaveLength(1);
    expect(arenas[0]).not.toBe(first);
    expect(arenas[0].pageCount).toBe(64);
  });
});
