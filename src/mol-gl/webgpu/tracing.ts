/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Camera } from '../../mol-canvas3d/camera';
import { TracingProps } from '../../mol-canvas3d/passes/tracing';
import { Mat4 } from '../../mol-math/linear-algebra';
import { RendererProps } from '../renderer';
import { GPUBufferUsage, GPUShaderStage, GPUTextureUsage } from './compat';
import { WebGPUContext } from './context';
import { createLightingData, WebGPULightingByteSize } from './lighting';
import { WebGPUTracingInput } from './tracing-input';
import { tracingShader } from './tracing-shader';

/** Screen-space diffuse path tracing with native PCG sampling and ping-pong accumulation. */
export class WebGPUTracing {
    private readonly settings: GPUBuffer;
    private readonly lights: GPUBuffer;
    private readonly layout: GPUBindGroupLayout;
    private readonly pipeline: GPUComputePipeline;
    private readonly targets: GPUTexture[] = [];
    private current = 0;
    private lastIteration = -1;
    private disposed = false;

    private constructor(private readonly context: WebGPUContext, module: GPUShaderModule) {
        const { device } = context;
        this.settings = device.createBuffer({ size: 208, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.lights = device.createBuffer({ size: WebGPULightingByteSize, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.layout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
            ...[2, 3, 4].map(binding => ({ binding, visibility: GPUShaderStage.COMPUTE, texture: { sampleType: 'float' as const } })),
            ...[5, 6].map(binding => ({ binding, visibility: GPUShaderStage.COMPUTE, texture: { sampleType: 'depth' as const } })),
            { binding: 7, visibility: GPUShaderStage.COMPUTE, texture: { sampleType: 'unfilterable-float' } },
            { binding: 8, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'write-only', format: 'rgba32float' } },
        ] });
        this.pipeline = device.createComputePipeline({ label: 'molstar-illumination-trace', layout: device.createPipelineLayout({ bindGroupLayouts: [this.layout] }), compute: { module, entryPoint: 'trace' } });
    }

    static async create(context: WebGPUContext) {
        const module = context.device.createShaderModule({ label: 'molstar-illumination-trace', code: tracingShader });
        const info = await module.getCompilationInfo();
        const errors = info.messages.filter(message => message.type === 'error');
        if (errors.length) throw new Error(`Native tracing shader compilation failed: ${errors.map(message => `${message.lineNum}:${message.linePos}: ${message.message}`).join('\n')}`);
        return new WebGPUTracing(context, module);
    }

    render(camera: Camera, renderer: RendererProps, props: TracingProps, input: WebGPUTracingInput['textures'], iteration: number, rays = props.rendersPerFrame[0], steps = props.steps, refineSteps = props.refineSteps, commandEncoder?: GPUCommandEncoder) {
        if (this.disposed) throw new Error('Native tracing has been disposed.');
        if (!Number.isInteger(iteration) || iteration < 0) throw new Error('Tracing iteration must be a nonnegative integer.');
        const { device } = this.context, { width, height } = input.shaded;
        if (this.targets[0]?.width !== width || this.targets[0]?.height !== height) {
            this.destroyTextures(); this.targets.length = 0;
            for (let i = 0; i < 2; i++) this.targets.push(device.createTexture({ label: 'molstar-illumination-accumulation', size: [width, height], format: 'rgba32float', usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC }));
            this.lastIteration = -1;
        }
        if (iteration !== 0 && iteration !== this.lastIteration + 1) throw new Error('Tracing iterations must be consecutive, starting at zero.');
        const data = new Float32Array(52), viewport = camera.viewport;
        data.set(camera.projection); data.set(Mat4.invert(Mat4(), camera.projection), 16);
        data.set([viewport.x, height - viewport.y - viewport.height, viewport.width, viewport.height], 32);
        data.set([width, height, camera.scale, iteration], 36);
        data.set([Math.max(1, Math.min(64, Math.round(rays))), Math.max(1, Math.min(1024, Math.round(steps))), Math.max(0, Math.min(8, Math.round(refineSteps))), Math.max(1, Math.min(32, Math.round(props.bounces)))], 40);
        data.set([props.rayDistance * camera.scale, props.minThickness * camera.scale, props.thicknessFactor, props.thickness * camera.scale], 44);
        data.set([props.thicknessMode === 'auto' ? 1 : 0, props.shadowEnable ? 1 : 0, props.shadowSoftness, props.shadowThickness * camera.scale], 48);
        device.queue.writeBuffer(this.settings, 0, data); device.queue.writeBuffer(this.lights, 0, createLightingData(camera, renderer, width, height));
        const previous = this.targets[this.current]; this.current = 1 - this.current;
        const output = this.targets[this.current];
        const bindings = device.createBindGroup({ layout: this.layout, entries: [
            { binding: 0, resource: { buffer: this.settings } }, { binding: 1, resource: { buffer: this.lights } },
            ...[input.shaded, input.normal, input.albedo, input.depth, input.backDepth, previous].map((texture, i) => ({ binding: i + 2, resource: texture.createView() })),
            { binding: 8, resource: output.createView() },
        ] });
        const encoder = commandEncoder ?? device.createCommandEncoder({ label: 'molstar-illumination-trace' });
        const pass = encoder.beginComputePass({ label: 'molstar-illumination-trace' });
        pass.setPipeline(this.pipeline); pass.setBindGroup(0, bindings); pass.dispatchWorkgroups(Math.ceil(width / 8), Math.ceil(height / 8)); pass.end();
        this.context.stats.computeDispatches++;
        if (!commandEncoder) device.queue.submit([encoder.finish()]);
        this.lastIteration = iteration;
        return output;
    }

    private destroyTextures() { for (const texture of this.targets) texture.destroy(); }
    dispose() { if (this.disposed) return; this.disposed = true; this.destroyTextures(); this.settings.destroy(); this.lights.destroy(); }
}
