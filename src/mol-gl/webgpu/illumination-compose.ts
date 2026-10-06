/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Camera } from '../../mol-canvas3d/camera';
import { IlluminationProps } from '../../mol-canvas3d/passes/illumination';
import { Mat4 } from '../../mol-math/linear-algebra';
import { Color } from '../../mol-util/color';
import { RendererProps } from '../renderer';
import { GPUBufferUsage, GPUShaderStage, GPUTextureUsage } from './compat';
import { WebGPUContext } from './context';
import { WebGPUTracingInput } from './tracing-input';

const shader = /* wgsl */ `
struct Settings { inverseProjection: mat4x4f, viewport: vec4f, dimensions: vec4f, background: vec4f, fog: vec4f };
@group(0) @binding(0) var<uniform> settings: Settings;
@group(0) @binding(1) var traced: texture_2d<f32>;
@group(0) @binding(2) var normals: texture_2d<f32>;
@group(0) @binding(3) var depth: texture_depth_2d;
fn coord(p: vec2i) -> vec2i { return clamp(p, vec2i(settings.viewport.xy), vec2i(settings.viewport.xy + settings.viewport.zw) - 1); }
fn denoise(p: vec2i) -> vec3f {
    let center = textureLoad(traced, p, 0);
    if (settings.fog.z <= 0.0) { return center.rgb; }
    let normal = textureLoad(normals, p, 0).xyz;
    let threshold = settings.fog.z; var sum = vec4f(0.0); var weightSum = 0.0;
    for (var x = -6; x <= 6; x++) { for (var y = -6; y <= 6; y++) {
        let d = vec2i(x, y); let q = coord(p + d);
        let n = textureLoad(normals, q, 0).xyz;
        let normalWeight = pow(clamp(dot(normal, n), 0.0, 1.0), 6.0);
        let sample = textureLoad(traced, q, 0); let delta = sample - center;
        let weight = exp(-f32(dot(d, d)) / 18.0) * (1.0 / (18.0 * 3.141592653589793)) * normalWeight
            * exp(-dot(delta, delta) / (2.0 * threshold * threshold)) * (0.3989422804014327 / threshold);
        sum += sample * weight; weightSum += weight;
    } }
    if (weightSum <= 0.000001) { return center.rgb; }
    return (sum / weightSum).rgb;
}
@vertex fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
    let vertices = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    return vec4f(vertices[i], 0.0, 1.0);
}
@fragment fn compose(@builtin(position) p: vec4f) -> @location(0) vec4f {
    let pixel = vec2i(p.xy);
    if (any(p.xy < settings.viewport.xy) || any(p.xy >= settings.viewport.xy + settings.viewport.zw) || textureLoad(depth, pixel, 0) == 1.0) {
        return settings.background;
    }
    let uv = (p.xy - settings.viewport.xy) / settings.viewport.zw;
    let view = settings.inverseProjection * vec4f(uv.x * 2.0 - 1.0, 1.0 - uv.y * 2.0, textureLoad(depth, pixel, 0) * 2.0 - 1.0, 1.0);
    var fog = 0.0;
    if (settings.fog.y > settings.fog.x) { fog = smoothstep(settings.fog.x, settings.fog.y, abs(view.z / view.w)); }
    let color = denoise(pixel);
    if (settings.dimensions.z > 0.0) { return vec4f(color * (1.0 - fog), 1.0 - fog); }
    return vec4f(mix(color, settings.background.rgb, fog), 1.0);
}
`;

/** Denoise and fog traced opaque lighting before normal transparent geometry and postprocessing. */
export class WebGPUIlluminationCompose {
    private readonly settings: GPUBuffer;
    private readonly layout: GPUBindGroupLayout;
    private readonly pipeline: GPURenderPipeline;
    private color?: GPUTexture;
    private disposed = false;

    private constructor(private readonly context: WebGPUContext, module: GPUShaderModule) {
        const { device } = context;
        this.settings = device.createBuffer({ size: 128, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.layout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'unfilterable-float' } },
            { binding: 2, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 3, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'depth' } },
        ] });
        this.pipeline = device.createRenderPipeline({ label: 'molstar-illumination-compose', layout: device.createPipelineLayout({ bindGroupLayouts: [this.layout] }),
            vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'compose', targets: [{ format: context.format }] },
        });
    }

    static async create(context: WebGPUContext) {
        const module = context.device.createShaderModule({ label: 'molstar-illumination-compose', code: shader });
        const info = await module.getCompilationInfo();
        const errors = info.messages.filter(message => message.type === 'error');
        if (errors.length) throw new Error(`Illumination composition shader failed: ${errors.map(message => message.message).join('\n')}`);
        return new WebGPUIlluminationCompose(context, module);
    }

    render(encoder: GPUCommandEncoder, traced: GPUTexture, input: WebGPUTracingInput['textures'], camera: Camera, props: IlluminationProps, renderer: RendererProps, transparentBackground: boolean, iteration: number, fullSampling = false, output?: GPUTexture) {
        if (this.disposed) throw new Error('Illumination composition has been disposed.');
        const { device } = this.context, { width, height } = traced;
        if (!output) {
            if (this.color?.width !== width || this.color?.height !== height) {
                this.color?.destroy(); this.color = device.createTexture({ label: 'molstar-illumination-color', size: [width, height], format: this.context.format, usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC });
            }
            output = this.color;
        }
        const data = new Float32Array(32), viewport = camera.viewport;
        data.set(Mat4.invert(Mat4(), camera.projection));
        data.set([viewport.x, height - viewport.y - viewport.height, viewport.width, viewport.height], 16);
        data.set([width, height, transparentBackground ? 1 : 0, 0], 20);
        const alpha = transparentBackground ? 0 : 1;
        data.set([...Color.toRgbNormalized(renderer.backgroundColor).map(v => v * alpha), alpha], 24);
        const progress = fullSampling ? 1 : Math.min(1, Math.max(0, iteration) / (Math.pow(2, props.maxIterations) / 2));
        const threshold = props.denoise ? props.denoiseThreshold[1] * (1 - progress) + props.denoiseThreshold[0] * progress : 0;
        data.set([camera.state.fog ? camera.fogNear : 0, camera.state.fog ? camera.fogFar : 0, threshold, 0], 28);
        device.queue.writeBuffer(this.settings, 0, data);
        const bindings = device.createBindGroup({ layout: this.layout, entries: [
            { binding: 0, resource: { buffer: this.settings } }, { binding: 1, resource: traced.createView() },
            { binding: 2, resource: input.normal.createView() }, { binding: 3, resource: input.depth.createView() },
        ] });
        const pass = encoder.beginRenderPass({ label: 'molstar-illumination-compose', colorAttachments: [{ view: output!.createView(), loadOp: 'clear', storeOp: 'store' }] });
        pass.setPipeline(this.pipeline); pass.setBindGroup(0, bindings); pass.draw(3); pass.end();
        return output!;
    }

    dispose() { if (this.disposed) return; this.disposed = true; this.color?.destroy(); this.settings.destroy(); }
}
