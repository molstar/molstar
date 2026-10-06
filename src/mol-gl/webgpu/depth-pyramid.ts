/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { WebGPUContext } from './context';
import { GPUShaderStage, GPUTextureUsage } from './compat';

const shader = /* wgsl */ `
@group(0) @binding(0) var opaque: texture_depth_2d;
@group(0) @binding(1) var transparent: texture_2d<f32>;
@group(0) @binding(2) var source: texture_2d<f32>;
@group(0) @binding(3) var outputDepth: texture_storage_2d<rgba32float, write>;
@compute @workgroup_size(8, 8) fn base(@builtin(global_invocation_id) id: vec3u) {
    let size = textureDimensions(outputDepth);
    if (any(id.xy >= size)) { return; }
    let p = min(vec2u((vec2f(id.xy) + 0.5) * vec2f(textureDimensions(opaque)) / vec2f(size)), textureDimensions(opaque) - 1u);
    let t = textureLoad(transparent, vec2i(p), 0).rg;
    textureStore(outputDepth, vec2i(id.xy), vec4f(textureLoad(opaque, vec2i(p), 0), t, 1.0));
}
@compute @workgroup_size(8, 8) fn reduce(@builtin(global_invocation_id) id: vec3u) {
    let size = textureDimensions(outputDepth);
    if (any(id.xy >= size)) { return; }
    // Match the nearest depth/alpha sampling used by the legacy half/quarter copy passes.
    let p = min(vec2u((vec2f(id.xy) + 0.5) * vec2f(textureDimensions(source)) / vec2f(size)), textureDimensions(source) - 1u);
    textureStore(outputDepth, vec2i(id.xy), textureLoad(source, vec2i(p), 0));
}
`;

/** Opaque depth and nearest transparent depth/alpha at SSAO, half and quarter resolution. */
export class WebGPUDepthPyramid {
    private texture?: GPUTexture;
    private readonly layout: GPUBindGroupLayout;
    private readonly base: GPUComputePipeline;
    private readonly reduce: GPUComputePipeline;
    constructor(private readonly context: WebGPUContext) {
        const { device } = context;
        this.layout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.COMPUTE, texture: { sampleType: 'depth' } },
            ...[1, 2].map(binding => ({ binding, visibility: GPUShaderStage.COMPUTE, texture: { sampleType: 'unfilterable-float' as const } })),
            { binding: 3, visibility: GPUShaderStage.COMPUTE, storageTexture: { access: 'write-only', format: 'rgba32float' } },
        ] });
        const layout = device.createPipelineLayout({ bindGroupLayouts: [this.layout] });
        const module = device.createShaderModule({ label: 'molstar-ssao-depth-pyramid', code: shader });
        this.base = device.createComputePipeline({ layout, compute: { module, entryPoint: 'base' } });
        this.reduce = device.createComputePipeline({ layout, compute: { module, entryPoint: 'reduce' } });
    }
    render(encoder: GPUCommandEncoder, opaque: GPUTexture, transparent: GPUTexture, width: number, height: number) {
        const levels = Math.min(3, Math.floor(Math.log2(Math.max(width, height))) + 1);
        if (!this.texture || this.texture.width !== width || this.texture.height !== height) {
            this.texture?.destroy();
            this.texture = this.context.device.createTexture({ label: 'molstar-ssao-depth-pyramid', size: { width, height }, mipLevelCount: levels, format: 'rgba32float', usage: GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC });
        }
        for (let level = 0; level < levels; level++) {
            const bindings = this.context.device.createBindGroup({ layout: this.layout, entries: [
                { binding: 0, resource: opaque.createView() }, { binding: 1, resource: transparent.createView() },
                { binding: 2, resource: level ? this.texture.createView({ baseMipLevel: level - 1, mipLevelCount: 1 }) : transparent.createView() },
                { binding: 3, resource: this.texture.createView({ baseMipLevel: level, mipLevelCount: 1 }) },
            ] });
            const pass = encoder.beginComputePass({ label: `molstar-ssao-depth-level-${level}` });
            pass.setPipeline(level ? this.reduce : this.base); pass.setBindGroup(0, bindings);
            pass.dispatchWorkgroups(Math.ceil(Math.max(1, width >> level) / 8), Math.ceil(Math.max(1, height >> level) / 8)); pass.end();
        }
        return this.texture;
    }
    dispose() { this.texture?.destroy(); this.texture = undefined; }
}
