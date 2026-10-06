/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { GPUBufferUsage, GPUMapMode } from './compat';

/** Read tightly packed rows from a native texture, including padded GPU copies. */
export async function readWebGPUTexture(device: GPUDevice, texture: GPUTexture, x: number, y: number, width: number, height: number, bytesPerPixel = 4, mipLevel = 0): Promise<Uint8Array> {
    if (![x, y, width, height, mipLevel].every(Number.isInteger) || mipLevel < 0 || mipLevel >= texture.mipLevelCount || x < 0 || y < 0 || width < 1 || height < 1 || x + width > Math.max(1, texture.width >> mipLevel) || y + height > Math.max(1, texture.height >> mipLevel)) {
        throw new Error('Readback rectangle is outside the WebGPU texture.');
    }
    const rowSize = width * bytesPerPixel;
    const bytesPerRow = Math.ceil(rowSize / 256) * 256;
    const buffer = device.createBuffer({ size: bytesPerRow * height, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
    try {
        const encoder = device.createCommandEncoder({ label: 'molstar-readback' });
        encoder.copyTextureToBuffer({ texture, mipLevel, origin: { x, y }, aspect: texture.format === 'depth32float' ? 'depth-only' : 'all' }, { buffer, bytesPerRow, rowsPerImage: height }, { width, height });
        device.queue.submit([encoder.finish()]);
        await buffer.mapAsync(GPUMapMode.READ);
        const mapped = new Uint8Array(buffer.getMappedRange());
        const result = new Uint8Array(rowSize * height);
        for (let row = 0; row < height; row++) result.set(mapped.subarray(row * bytesPerRow, row * bytesPerRow + rowSize), row * rowSize);
        buffer.unmap();
        return result;
    } finally {
        buffer.destroy();
    }
}
