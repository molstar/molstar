/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Camera } from '../../mol-canvas3d/camera';
import { WebGPUContext } from './context';
import { GPUShaderStage, GPUTextureUsage } from './compat';

export const weightedTransparencyShader = /* wgsl */ `
struct WeightedFragment { @location(0) color: vec4f, @location(1) weight: vec4f, @location(2) emission: vec4f };
fn weightedFragment(color: vec4f, emission: vec4f, depth: f32) -> WeightedFragment {
    let weight = color.a * clamp((1.0 - depth) * (1.0 - depth), 0.01, 1.0);
    return WeightedFragment(vec4f(color.rgb * weight, color.a),
        vec4f(color.a * weight, emission.a * weight, 0.0, emission.a), vec4f(emission.rgb * weight, 0.0));
}
`;

const shader = /* wgsl */ `
@group(0) @binding(0) var accum: texture_2d<f32>;
@group(0) @binding(1) var weights: texture_2d<f32>;
@group(0) @binding(2) var emission: texture_2d<f32>;
@vertex fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
    let p = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    return vec4f(p[i], 0.0, 1.0);
}
struct Output { @location(0) color: vec4f, @location(1) emission: vec4f };
@fragment fn fs(@builtin(position) p: vec4f) -> Output {
    let c = textureLoad(accum, vec2i(p.xy), 0); let w = textureLoad(weights, vec2i(p.xy), 0);
    let e = textureLoad(emission, vec2i(p.xy), 0);
    let alpha = 1.0 - c.a; let coverage = 1.0 - w.a;
    return Output(vec4f(c.rgb / max(w.r, 0.00000001) * alpha, alpha), vec4f(e.rgb / max(w.g, 0.00000001) * coverage, coverage));
}
`;

/** Portable weighted color/emission accumulation; opaque and selection depths remain independent. */
export class WebGPUWeightedTransparency {
    private textures: GPUTexture[] = [];
    private readonly layout: GPUBindGroupLayout;
    private readonly pipeline: GPURenderPipeline;
    private readonly colorPipeline: GPURenderPipeline;
    private readonly noEmissionPipeline: GPURenderPipeline;
    static targets(emission = true): GPUColorTargetState[] {
        const blend: GPUBlendState = { color: { srcFactor: 'one', dstFactor: 'one' }, alpha: { srcFactor: 'zero', dstFactor: 'one-minus-src-alpha' } };
        return [{ format: 'rgba16float', blend }, { format: 'rgba16float', blend },
            { format: 'rgba16float', writeMask: emission ? 15 : 0, blend: { color: { srcFactor: 'one', dstFactor: 'one' }, alpha: { srcFactor: 'one', dstFactor: 'one' } } }];
    }
    constructor(private readonly context: WebGPUContext) {
        const { device } = context;
        this.layout = device.createBindGroupLayout({ entries: [0, 1, 2].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' as const } })) });
        const module = device.createShaderModule({ label: 'molstar-weighted-transparency', code: shader });
        const blend: GPUBlendState = { color: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' }, alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' } };
        const layout = device.createPipelineLayout({ bindGroupLayouts: [this.layout] });
        const create = (mask: number) => device.createRenderPipeline({ layout,
            vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'fs', targets: [{ format: context.format, blend }, { format: 'rgba16float', blend, writeMask: mask }] } });
        this.pipeline = create(15); this.colorPipeline = create(0); this.noEmissionPipeline = create(8);
    }
    begin(encoder: GPUCommandEncoder, depth: GPUTexture, camera: Camera) {
        if (this.textures[0]?.width !== depth.width || this.textures[0]?.height !== depth.height) {
            this.destroyTextures();
            this.textures = [0, 1, 2].map(i => this.context.device.createTexture({ label: `molstar-weighted-${i}`, size: [depth.width, depth.height], format: 'rgba16float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING }));
        }
        const pass = encoder.beginRenderPass({ label: 'molstar-weighted-accumulation',
            colorAttachments: this.textures.map((texture, i) => ({ view: texture.createView(), loadOp: 'clear' as const, storeOp: 'store' as const, clearValue: { r: 0, g: 0, b: 0, a: i < 2 ? 1 : 0 } })),
            depthStencilAttachment: { view: depth.createView(), depthReadOnly: true },
        });
        const v = camera.viewport; pass.setViewport(v.x, depth.height - v.y - v.height, v.width, v.height, 0, 1);
        return pass;
    }
    resolve(encoder: GPUCommandEncoder, color: GPUTexture, emission: GPUTexture, camera: Camera, transparentColor?: GPUTexture, emitTransparent = true) {
        const group = this.context.device.createBindGroup({ layout: this.layout, entries: this.textures.map((texture, binding) => ({ binding, resource: texture.createView() })) });
        const run = (output: GPUTexture, clear: boolean) => {
            const pass = encoder.beginRenderPass({ label: 'molstar-weighted-resolve', colorAttachments: [output, emission].map((texture, i) => ({ view: texture.createView(), loadOp: clear && i === 0 ? 'clear' as const : 'load' as const, storeOp: 'store' as const })) });
            const v = camera.viewport; pass.setViewport(v.x, output.height - v.y - v.height, v.width, v.height, 0, 1);
            pass.setPipeline(clear ? this.colorPipeline : emitTransparent ? this.pipeline : this.noEmissionPipeline); pass.setBindGroup(0, group); pass.draw(3); pass.end();
        };
        // Render the transparent-only color first; the final resolve supplies emission once.
        if (transparentColor) run(transparentColor, true);
        run(color, false);
    }
    private destroyTextures() { for (const texture of this.textures) texture.destroy(); this.textures = []; }
    dispose() { this.destroyTextures(); }
}
