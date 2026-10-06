/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 */

import { Subject } from 'rxjs';
import { readWebGPUTexture } from './readback';
import { GPUTextureUsage } from './compat';

/** Native WebGPU device and canvas ownership. Does not acquire a WebGL context. */
export class WebGPUContext {
    readonly lost = new Subject<GPUDeviceLostInfo>();
    readonly errors = new Subject<GPUError>();
    readonly stats = { computeDispatches: 0, marchingCubesDispatches: 0 };
    readonly format: GPUTextureFormat;
    private disposed = false;

    private constructor(readonly canvas: { width: number, height: number }, readonly adapter: GPUAdapter, readonly device: GPUDevice, readonly context: GPUCanvasContext | undefined, gpu: GPU) {
        this.format = gpu.getPreferredCanvasFormat();
        context?.configure({ device, format: this.format, alphaMode: 'premultiplied', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_DST });
        device.addEventListener('uncapturederror', this.handleError);
        device.lost.then(info => {
            if (!this.disposed) this.lost.next(info);
        });
    }

    private handleError = (event: GPUUncapturedErrorEvent) => this.errors.next(event.error);

    static async create(canvas: HTMLCanvasElement, options: GPURequestAdapterOptions = {}, gpu: GPU | undefined = typeof navigator === 'undefined' ? undefined : navigator.gpu): Promise<WebGPUContext> {
        if (!gpu) throw new Error('WebGPU is unavailable. Use a WebGPU browser on HTTPS or localhost.');
        const adapter = await gpu.requestAdapter({ powerPreference: 'high-performance', ...options });
        if (!adapter) throw new Error('No WebGPU adapter is available.');
        const device = await adapter.requestDevice();
        try {
            const context = canvas.getContext('webgpu') as unknown as GPUCanvasContext | null;
            if (!context) throw new Error('Could not create a WebGPU canvas context.');
            return new WebGPUContext(canvas, adapter, device, context, gpu);
        } catch (error) {
            device.destroy();
            throw error;
        }
    }

    /** Offscreen native textures for Node.js or other environments without a DOM. */
    static async createHeadless(size: { width: number, height: number }, gpu: GPU, options: GPURequestAdapterOptions = {}) {
        if (![size.width, size.height].every(v => Number.isInteger(v) && v > 0)) throw new Error('Headless WebGPU dimensions must be positive integers.');
        const adapter = await gpu.requestAdapter({ powerPreference: 'high-performance', ...options });
        if (!adapter) throw new Error('No headless WebGPU adapter is available.');
        const device = await adapter.requestDevice();
        if (Math.max(size.width, size.height) > device.limits.maxTextureDimension2D) { device.destroy(); throw new Error('Headless dimensions exceed the WebGPU texture limit.'); }
        return new WebGPUContext({ ...size }, adapter, device, undefined, gpu);
    }

    /** Allocate a GPU buffer with a padded size, preserving typed-array offsets. */
    createBuffer(data: ArrayBufferView, usage: GPUBufferUsageFlags, label?: string): GPUBuffer {
        if (this.disposed) throw new Error('WebGPU context has been disposed.');
        const size = Math.max(4, Math.ceil(data.byteLength / 4) * 4);
        const buffer = this.device.createBuffer({ size, usage, label, mappedAtCreation: true });
        new Uint8Array(buffer.getMappedRange()).set(new Uint8Array(data.buffer, data.byteOffset, data.byteLength));
        buffer.unmap();
        return buffer;
    }

    /** WebGPU rows are aligned to 256 bytes; callers receive tightly packed rows. */
    async readTexture(texture: GPUTexture, x: number, y: number, width: number, height: number, bytesPerPixel = 4, mipLevel = 0): Promise<Uint8Array> {
        if (this.disposed) throw new Error('WebGPU context has been disposed.');
        return readWebGPUTexture(this.device, texture, x, y, width, height, bytesPerPixel, mipLevel);
    }

    dispose() {
        if (this.disposed) return;
        this.disposed = true;
        this.device.removeEventListener('uncapturederror', this.handleError);
        this.context?.unconfigure();
        this.device.destroy();
        this.lost.complete();
        this.errors.complete();
    }
}
