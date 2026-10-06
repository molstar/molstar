/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Camera } from '../../mol-canvas3d/camera';
import { BloomProps } from '../../mol-canvas3d/passes/bloom';
import { GPUBufferUsage, GPUShaderStage, GPUTextureUsage } from './compat';
import { WebGPUContext } from './context';

export const bloomShader = /* wgsl */ `
struct Settings { dimensions: vec4f, config: vec4f, direction: vec4f, viewport: vec4f };
@group(0) @binding(0) var<uniform> settings: Settings;
@group(0) @binding(1) var color: texture_2d<f32>;
@group(0) @binding(2) var emissive: texture_2d<f32>;
@group(0) @binding(3) var depth: texture_depth_2d;
@group(0) @binding(4) var picking: texture_2d<u32>;
@group(0) @binding(5) var linearSampler: sampler;
@group(0) @binding(6) var mip1: texture_2d<f32>;
@group(0) @binding(7) var mip2: texture_2d<f32>;
@group(0) @binding(8) var mip3: texture_2d<f32>;
@group(0) @binding(9) var mip4: texture_2d<f32>;
@group(0) @binding(10) var mip5: texture_2d<f32>;
@vertex fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
    let positions = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    return vec4f(positions[i], 0.0, 1.0);
}
fn outside(uv: vec2f) -> bool { return any(uv < settings.viewport.xy) || any(uv >= settings.viewport.xy + settings.viewport.zw); }
fn sampleInput(uv: vec2f) -> vec4f {
    let halfPixel = 0.5 / vec2f(textureDimensions(color));
    let center = settings.viewport.xy + settings.viewport.zw * 0.5;
    let low = min(center, settings.viewport.xy + halfPixel);
    let high = max(center, settings.viewport.xy + settings.viewport.zw - halfPixel);
    return textureSampleLevel(color, linearSampler, clamp(uv, low, high), 0.0);
}
@fragment fn seed(@builtin(position) p: vec4f) -> @location(0) vec4f {
    let uv = p.xy / settings.dimensions.xy; if (outside(uv)) { return vec4f(0.0); }
    let pixel = vec2i(p.xy); var d = textureLoad(depth, pixel, 0); let pick = textureLoad(picking, pixel, 0);
    if (pick.x > 0u) { d = min(d, bitcast<f32>(pick.w)); }
    if (d >= 0.99999994) { return vec4f(0.0); }
    if (settings.dimensions.w > 0.0) { return textureLoad(emissive, pixel, 0); }
    let texel = textureLoad(color, pixel, 0); let luminosity = dot(texel.rgb, vec3f(0.299, 0.587, 0.114));
    return texel * smoothstep(settings.dimensions.z, settings.dimensions.z + 0.01, luminosity);
}
@fragment fn blur(@builtin(position) p: vec4f) -> @location(0) vec4f {
    let uv = p.xy / settings.dimensions.xy; if (outside(uv)) { return vec4f(0.0); }
    var sum = sampleInput(uv); var weights = 1.0;
    let radius = settings.config.z;
    for (var i = 1u; i < u32(radius); i++) {
        let x = f32(i); let weight = exp(-0.5 * x * x / (radius * radius));
        let offset = settings.direction.xy / settings.dimensions.xy * x;
        sum += (sampleInput(uv + offset) + sampleInput(uv - offset)) * weight; weights += 2.0 * weight;
    }
    return sum / weights;
}
fn factor(value: f32) -> f32 { return mix(value, 1.2 - value, settings.config.y); }
@fragment fn composite(@builtin(position) p: vec4f) -> @location(0) vec4f {
    let uv = p.xy / settings.dimensions.xy; if (outside(uv)) { return vec4f(0.0); }
    return settings.config.x * (factor(1.0) * textureSampleLevel(mip1, linearSampler, uv, 0.0) +
        factor(0.8) * textureSampleLevel(mip2, linearSampler, uv, 0.0) + factor(0.6) * textureSampleLevel(mip3, linearSampler, uv, 0.0) +
        factor(0.4) * textureSampleLevel(mip4, linearSampler, uv, 0.0) + factor(0.2) * textureSampleLevel(mip5, linearSampler, uv, 0.0));
}
`;

/** Five progressively smaller Gaussian blur levels, matching Mol*'s bloom controls. */
export class WebGPUBloom {
    private readonly layout: GPUBindGroupLayout;
    private readonly pipelines: Record<'seed' | 'blur' | 'composite', GPURenderPipeline>;
    private readonly sampler: GPUSampler;
    private readonly uniforms = new Map<string, GPUBuffer>();
    private readonly horizontal: GPUTexture[] = [];
    private readonly vertical: GPUTexture[] = [];
    private seed?: GPUTexture;
    private output?: GPUTexture;

    constructor(private readonly context: WebGPUContext, shader: GPUShaderModule) {
        const { device } = context;
        this.layout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            ...[1, 2, 6, 7, 8, 9, 10].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' as const } })),
            { binding: 3, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'depth' } },
            { binding: 4, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'uint' } },
            { binding: 5, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
        ] });
        const layout = device.createPipelineLayout({ bindGroupLayouts: [this.layout] });
        const pipeline = (entryPoint: string) => device.createRenderPipeline({ label: `molstar-bloom-${entryPoint}`, layout,
            vertex: { module: shader, entryPoint: 'vs' }, fragment: { module: shader, entryPoint, targets: [{ format: 'rgba16float' }] }, primitive: { topology: 'triangle-list' } });
        this.pipelines = { seed: pipeline('seed'), blur: pipeline('blur'), composite: pipeline('composite') };
        this.sampler = device.createSampler({ minFilter: 'linear', magFilter: 'linear' });
    }

    private resize(width: number, height: number) {
        if (this.output?.width === width && this.output.height === height) return;
        this.destroyTextures();
        const texture = (width: number, height: number) => this.context.device.createTexture({ size: { width, height }, format: 'rgba16float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING });
        this.seed = texture(width, height); this.output = texture(width, height);
        for (let i = 0; i < 5; i++) {
            width = Math.max(1, Math.round(width / 2)); height = Math.max(1, Math.round(height / 2));
            this.horizontal.push(texture(width, height)); this.vertical.push(texture(width, height));
        }
    }

    render(encoder: GPUCommandEncoder, color: GPUTexture, emissive: GPUTexture, depth: GPUTexture, picking: GPUTexture, camera: Camera, props: BloomProps) {
        this.resize(color.width, color.height);
        const run = (name: string, pipeline: GPURenderPipeline, input: GPUTexture, target: GPUTexture, radius = 0, direction: [number, number] = [0, 0], mips?: GPUTexture[]) => {
            let uniform = this.uniforms.get(name);
            if (!uniform) { uniform = this.context.device.createBuffer({ size: 64, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST }); this.uniforms.set(name, uniform); }
            const v = camera.viewport;
            const data = new Float32Array([target.width, target.height, props.threshold, props.mode === 'emissive' ? 1 : 0,
                props.strength, props.radius, radius, 0, ...direction, 0, 0,
                v.x / color.width, (color.height - v.y - v.height) / color.height, v.width / color.width, v.height / color.height]);
            this.context.device.queue.writeBuffer(uniform, 0, data);
            const bindings = this.context.device.createBindGroup({ layout: this.layout, entries: [
                { binding: 0, resource: { buffer: uniform } }, { binding: 1, resource: input.createView() },
                { binding: 2, resource: emissive.createView() }, { binding: 3, resource: depth.createView() }, { binding: 4, resource: picking.createView() }, { binding: 5, resource: this.sampler },
                ...[6, 7, 8, 9, 10].map((binding, i) => ({ binding, resource: (mips?.[i] ?? color).createView() })),
            ] });
            const pass = encoder.beginRenderPass({ label: `molstar-${name}`, colorAttachments: [{ view: target.createView(), loadOp: 'clear', storeOp: 'store' }] });
            pass.setPipeline(pipeline); pass.setBindGroup(0, bindings); pass.draw(3); pass.end();
        };
        run('bloom-seed', this.pipelines.seed, color, this.seed!);
        for (let i = 0; i < 5; i++) {
            run(`bloom-horizontal-${i}`, this.pipelines.blur, i === 0 ? this.seed! : this.vertical[i - 1], this.horizontal[i], 3 + i * 2, [1, 0]);
            run(`bloom-vertical-${i}`, this.pipelines.blur, this.horizontal[i], this.vertical[i], 3 + i * 2, [0, 1]);
        }
        run('bloom-composite', this.pipelines.composite, color, this.output!, 0, [0, 0], this.vertical);
        return this.output!;
    }

    private destroyTextures() { this.seed?.destroy(); this.output?.destroy(); for (const t of [...this.horizontal, ...this.vertical]) t.destroy(); this.horizontal.length = 0; this.vertical.length = 0; }
    dispose() { this.destroyTextures(); for (const u of this.uniforms.values()) u.destroy(); this.uniforms.clear(); }
}
