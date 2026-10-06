/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Camera } from '../../mol-canvas3d/camera';
import { JitterVectors, MultiSampleProps } from '../../mol-canvas3d/passes/multi-sample';
import { WebGPUContext } from './context';
import { GPUBufferUsage, GPUShaderStage, GPUTextureUsage } from './compat';

const shader = /* wgsl */ `
@group(0) @binding(0) var<uniform> weight: vec4f;
@group(0) @binding(1) var source: texture_2d<f32>;
@group(0) @binding(2) var hold: texture_2d<f32>;
@vertex fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
    let vertices = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    return vec4f(vertices[i], 0.0, 1.0);
}
@fragment fn accumulate(@builtin(position) p: vec4f) -> @location(0) vec4f { return textureLoad(source, vec2i(p.xy), 0) * weight.x; }
@fragment fn compose(@builtin(position) p: vec4f) -> @location(0) vec4f { return textureLoad(source, vec2i(p.xy), 0) + textureLoad(hold, vec2i(p.xy), 0) * weight.y; }
`;

/** Supersampling with the canonical subpixel jitter vectors and premultiplied fp16 accumulation. */
export class WebGPUMultiSample {
    private readonly layout: GPUBindGroupLayout;
    private readonly uniform: GPUBuffer;
    private readonly accumulate: GPURenderPipeline;
    private readonly compose: GPURenderPipeline;
    private accumulation?: GPUTexture;
    private hold?: GPUTexture;
    private color?: GPUTexture;
    private picking?: GPUTexture;
    private progress = -1;
    private key = '';
    get needsFrame() { return this.progress >= 0; }
    reset() { this.progress = -1; this.key = ''; }
    constructor(private readonly context: WebGPUContext) {
        const { device } = context;
        this.layout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            ...[1, 2].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' as const } })),
        ] });
        this.uniform = device.createBuffer({ size: 16, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        const module = device.createShaderModule({ label: 'molstar-multi-sample', code: shader });
        const layout = device.createPipelineLayout({ bindGroupLayouts: [this.layout] });
        this.accumulate = device.createRenderPipeline({ layout, vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'accumulate', targets: [{ format: 'rgba16float', blend: { color: { srcFactor: 'one', dstFactor: 'one' }, alpha: { srcFactor: 'one', dstFactor: 'one' } } }] } });
        this.compose = device.createRenderPipeline({ layout, vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'compose', targets: [{ format: context.format }] } });
    }
    private setSize(width: number, height: number) {
        if (this.color?.width === width && this.color.height === height) return;
        this.destroyTextures(); this.reset();
        const usage = GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC | GPUTextureUsage.COPY_DST;
        this.accumulation = this.context.device.createTexture({ size: [width, height], format: 'rgba16float', usage });
        this.hold = this.context.device.createTexture({ size: [width, height], format: this.context.format, usage });
        this.color = this.context.device.createTexture({ size: [width, height], format: this.context.format, usage });
        this.picking = this.context.device.createTexture({ size: [width, height], format: 'rgba32uint', usage: GPUTextureUsage.COPY_DST | GPUTextureUsage.COPY_SRC });
    }
    private copy(source: GPUTexture, destination: GPUTexture) {
        const encoder = this.context.device.createCommandEncoder(); encoder.copyTextureToTexture({ texture: source }, { texture: destination }, [source.width, source.height]); this.context.device.queue.submit([encoder.finish()]);
    }
    private pass(source: GPUTexture, output: GPUTexture, pipeline: GPURenderPipeline, first: boolean, sampleWeight: number, holdWeight = 0) {
        const { device } = this.context;
        device.queue.writeBuffer(this.uniform, 0, new Float32Array([sampleWeight, holdWeight, 0, 0]));
        const bindings = device.createBindGroup({ layout: this.layout, entries: [{ binding: 0, resource: { buffer: this.uniform } }, { binding: 1, resource: source.createView() }, { binding: 2, resource: this.hold!.createView() }] });
        const encoder = device.createCommandEncoder({ label: 'molstar-multi-sample' });
        const pass = encoder.beginRenderPass({ colorAttachments: [{ view: output.createView(), loadOp: first ? 'clear' : 'load', storeOp: 'store' }] });
        pass.setPipeline(pipeline); pass.setBindGroup(0, bindings); pass.draw(3); pass.end(); device.queue.submit([encoder.finish()]);
    }
    render(camera: Camera, props: MultiSampleProps, changed: boolean, forceOn: boolean, width: number, height: number, picking: GPUTexture, draw: (offset?: readonly number[]) => GPUTexture) {
        this.setSize(width, height);
        const level = Math.max(0, Math.min(5, Math.round(props.sampleLevel))), offsets = JitterVectors[level];
        const key = JSON.stringify([props, camera.viewport]);
        changed ||= key !== this.key; this.key = key;
        const full = props.mode === 'on' || forceOn;
        const original = { ...camera.viewOffset };
        try {
            if (changed || full || !this.color) {
                const baseline = draw(); this.copy(baseline, this.hold!); this.copy(picking, this.picking!);
                this.progress = 0;
                if (!full) { this.copy(baseline, this.color!); return this.color!; }
            }
            if (this.progress < 0) return this.color!;
            const end = Math.min(offsets.length, this.progress + (full ? offsets.length : Math.pow(2, Math.max(0, level - 2))));
            for (; this.progress < end; this.progress++) {
                const offset = offsets[this.progress];
                camera.viewOffset.enabled = true;
                if (original.enabled) {
                    Camera.setViewOffset(camera.viewOffset, original.fullWidth, original.fullHeight,
                        original.offsetX + offset[0] * original.width / camera.viewport.width,
                        original.offsetY + offset[1] * original.height / camera.viewport.height, original.width, original.height);
                } else {
                    Camera.setViewOffset(camera.viewOffset, camera.viewport.width, camera.viewport.height, offset[0], offset[1], camera.viewport.width, camera.viewport.height);
                }
                camera.update();
                const source = draw(props.reuseOcclusion ? offset : undefined);
                const sampleWeight = 1 / offsets.length + (full ? (1 / 32) * (-0.5 + (this.progress + 0.5) / offsets.length) : 0);
                this.pass(source, this.accumulation!, this.accumulate, this.progress === 0, sampleWeight);
            }
            this.pass(this.accumulation!, this.color!, this.compose, true, 1, 1 - this.progress / offsets.length);
            if (this.progress === offsets.length) this.progress = -1;
            return this.color!;
        } finally {
            Object.assign(camera.viewOffset, original); camera.update();
            this.copy(this.picking!, picking);
        }
    }
    /** Progressive illumination keeps the canonical first frame while averaging jittered iterations. */
    renderIllumination(camera: Camera, props: MultiSampleProps, iteration: number, maxIterations: number, width: number, height: number, picking: GPUTexture, draw: (refresh: boolean) => GPUTexture) {
        this.setSize(width, height);
        if (iteration === 0) {
            const baseline = draw(true);
            this.copy(baseline, this.hold!); this.copy(baseline, this.color!); this.copy(picking, this.picking!);
            return this.color!;
        }
        const offsets = JitterVectors[Math.max(0, Math.min(5, Math.round(props.sampleLevel)))];
        const index = Math.min(offsets.length - 1, Math.floor(iteration * offsets.length / maxIterations));
        const previousIndex = Math.floor((iteration - 1) * offsets.length / maxIterations);
        const offset = offsets[index], original = { ...camera.viewOffset };
        try {
            camera.viewOffset.enabled = true;
            if (original.enabled) {
                Camera.setViewOffset(camera.viewOffset, original.fullWidth, original.fullHeight,
                    original.offsetX + offset[0] * original.width / camera.viewport.width,
                    original.offsetY + offset[1] * original.height / camera.viewport.height, original.width, original.height);
            } else {
                Camera.setViewOffset(camera.viewOffset, camera.viewport.width, camera.viewport.height, offset[0], offset[1], camera.viewport.width, camera.viewport.height);
            }
            camera.update();
            const source = draw(iteration === 1 || index !== previousIndex);
            this.pass(source, this.accumulation!, this.accumulate, iteration === 1, 1 / maxIterations);
            this.pass(this.accumulation!, this.color!, this.compose, true, 1, 1 - iteration / maxIterations);
            return this.color!;
        } finally {
            Object.assign(camera.viewOffset, original); camera.update(); this.copy(this.picking!, picking);
        }
    }

    private destroyTextures() { for (const texture of [this.accumulation, this.hold, this.color, this.picking]) texture?.destroy(); }
    dispose() { this.destroyTextures(); this.uniform.destroy(); }
}
