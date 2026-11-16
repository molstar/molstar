/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { createWebGPUContext } from '../webgpu/context.js';
import { createStubDevice, stubCalls } from './stub-gpu.js';

const Size = 48;

/** A WebGPU context on the stub device and four rgba8 textures, as a pick pass reads them */
function createCase() {
  const device = createStubDevice();
  // an own method of the stub (`stubCalls` counts the proxied ones)
  const createBuffer = vi.spyOn(device, 'createBuffer');
  const ctx = createWebGPUContext(device, { width: Size, height: Size });
  const textures = Array.from({ length: 4 }, () => {
    const t = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    t.define(Size, Size);
    return t;
  });
  const arrays = textures.map(() => new Uint8Array(Size * Size * 4));
  const buffers = () => createBuffer.mock.calls.length;
  return { ctx, textures, arrays, buffers };
}

const count = (name: string) => stubCalls.get(name) ?? 0;

describe('a WebGPU texture reader', () => {
  it('reads the textures of a read with one staging buffer, command buffer and submit, the buffer kept', async () => {
    const { ctx, textures, arrays, buffers: madeBuffers } = createCase();
    const reader = ctx.createTextureReader();
    const read = async () => {
      reader.start(textures, 0, 0, 30, 20);
      for (let k = 0; k < 10 && reader.check(arrays) === 'pending'; ++k) await Promise.resolve();
    };
    await read();
    const buffers = madeBuffers();
    const encoders = count('createCommandEncoder');
    const copies = count('copyTextureToBuffer');
    const reads = 8;
    for (let i = 0; i < reads; ++i) await read();
    // no staging buffer made again; one command buffer of the four copies a read
    expect(madeBuffers() - buffers).toBe(0);
    expect(count('createCommandEncoder') - encoders).toBe(reads);
    expect(count('copyTextureToBuffer') - copies).toBe(4 * reads);
    // the staging buffer counted while kept, released with the reader
    expect(reader.getByteCount()).toBeGreaterThan(0);
    expect(ctx.stats.resourceCounts.pixelPack).toBe(1);
    reader.dispose();
    expect(ctx.stats.resourceCounts.pixelPack).toBe(0);
    ctx.destroy();
  });

  it('drops a read started before the last one, one staging buffer kept for the next', async () => {
    const { ctx, textures, arrays, buffers: madeBuffers } = createCase();
    const reader = ctx.createTextureReader();
    const before = madeBuffers();
    reader.start(textures, 0, 0, 16, 16);
    reader.start(textures, 0, 0, 16, 16);
    for (let k = 0; k < 10 && reader.check(arrays) === 'pending'; ++k) await Promise.resolve();
    // two reads in flight at once made two buffers, one kept
    expect(madeBuffers() - before).toBe(2);
    expect(ctx.stats.resourceCounts.pixelPack).toBe(1);
    const buffers = madeBuffers();
    reader.start(textures, 0, 0, 16, 16);
    expect(madeBuffers() - buffers).toBe(0);
    reader.dispose();
    for (let k = 0; k < 10; ++k) await Promise.resolve();
    expect(ctx.stats.resourceCounts.pixelPack).toBe(0);
    ctx.destroy();
  });
});
