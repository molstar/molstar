/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { TextureImage, TextureVolume } from '../renderable/util';
import type { Texture } from '../webgl/texture';
import { GPUTextureUsage } from './compat';
import { readWebGPUTexture } from './readback';
import { idFactory } from '../../mol-util/id-factory';

const nextId = idFactory();
export interface VolumeTextureData {
    array: Uint8Array | Uint16Array | Float32Array
    width: number
    height: number
    depth: number
}

/** RGBA data or a native GPU texture. No GL resource is allocated. */
export class WebGPUTextureData implements Texture {
    static disposeTextures(meta: Record<string, unknown>) {
        for (const value of Object.values(meta)) if (value instanceof WebGPUTextureData) value.destroy();
    }
    readonly id = nextId();
    readonly target = 0;
    readonly format = 0;
    readonly internalFormat = 0;
    readonly type = 0;
    readonly filter = 0;
    private image?: VolumeTextureData;
    private destroyed = false;
    private revision = 0;
    get version() { return this.revision; }
    private gpuImage?: { device: GPUDevice, texture: GPUTexture };
    get native() { return this.destroyed ? undefined : this.gpuImage; }
    get data() {
        if (this.gpuImage && !this.image?.array.length) throw new Error('Native GPU texture data requires explicit readback.');
        if (!this.image || this.destroyed) throw new Error('WebGPU texture data is unavailable.');
        return this.image;
    }
    /** Export CPU-owned data or explicitly read a native RGBA8 grid. */
    async readData(): Promise<VolumeTextureData> {
        if (!this.gpuImage) return this.data;
        if (this.destroyed) throw new Error('WebGPU texture data is unavailable.');
        const { device, texture } = this.gpuImage;
        if (!(texture.usage & GPUTextureUsage.COPY_SRC)) throw new Error('Native texture readback requires COPY_SRC usage.');
        const bytes = await readWebGPUTexture(device, texture, 0, 0, texture.width, texture.height, texture.format === 'rgba32float' ? 16 : 4);
        const array = texture.format === 'rgba32float' ? new Float32Array(bytes.buffer) : bytes;
        return { array, width: texture.width, height: texture.height, depth: 1 };
    }
    getWidth() { return this.image?.width ?? 0; }
    getHeight() { return this.image?.height ?? 0; }
    getDepth() { return this.image?.depth ?? 0; }
    getByteCount() { return this.gpuImage ? this.getWidth() * this.getHeight() * this.getDepth() * (this.gpuImage.texture.format === 'rgba32float' ? 16 : 4) : this.image?.array.byteLength ?? 0; }
    define(_width: number, _height: number, _depth?: number) { throw new Error('Upload texture data with load before WebGPU rendering.'); }
    load(image: TextureImage<any> | TextureVolume<any> | HTMLImageElement, sub = false) {
        if (this.destroyed) throw new Error('Cannot upload a destroyed WebGPU texture.');
        if (sub || !('array' in image)) throw new Error('WebGPU texture data requires a complete array.');
        const data: VolumeTextureData = 'depth' in image ? image : { ...image, depth: 1 };
        if (![data.width, data.height, data.depth].every(n => Number.isInteger(n) && n > 0) || data.array.length < data.width * data.height * data.depth * 4) throw new Error('Invalid WebGPU texture dimensions or data length.');
        this.gpuImage?.texture.destroy(); this.gpuImage = undefined;
        this.image = data; this.revision++;
    }
    /** Take ownership of a texture produced by native compute. */
    loadGPU(device: GPUDevice, texture: GPUTexture, mirror?: VolumeTextureData) {
        if (this.destroyed) throw new Error('Cannot upload a destroyed WebGPU texture.');
        if (texture.dimension !== '2d' || texture.depthOrArrayLayers !== 1 || !['rgba8unorm', 'rgba32float'].includes(texture.format) || !(texture.usage & GPUTextureUsage.TEXTURE_BINDING)) throw new Error('Native textures require 2D RGBA8 or RGBA32F data.');
        if (mirror && (mirror.width !== texture.width || mirror.height !== texture.height || mirror.depth !== 1 || mirror.array.length !== texture.width * texture.height * 4 || (texture.format === 'rgba32float' ? !(mirror.array instanceof Float32Array) : !(mirror.array instanceof Uint8Array)))) throw new Error('Native texture mirror must match the GPU texture.');
        if (this.gpuImage?.texture !== texture) this.gpuImage?.texture.destroy();
        this.gpuImage = { device, texture }; this.revision++;
        this.image = mirror ?? { width: texture.width, height: texture.height, depth: 1, array: new Uint8Array(0) };
    }
    mipmap() { throw new Error('WebGPU volume data does not expose WebGL mipmaps.'); }
    bind() { throw new Error('WebGPU volume data cannot be bound to WebGL.'); }
    unbind() { throw new Error('WebGPU volume data cannot be bound to WebGL.'); }
    attachFramebuffer() { throw new Error('WebGPU volume data cannot be attached to a WebGL framebuffer.'); }
    detachFramebuffer() { throw new Error('WebGPU volume data cannot be attached to a WebGL framebuffer.'); }
    reset() { }
    destroy() { if (this.destroyed) return; this.revision++; this.destroyed = true; this.image = undefined; this.gpuImage?.texture.destroy(); this.gpuImage = undefined; }
}
