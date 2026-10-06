/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Camera } from '../../mol-canvas3d/camera';
import { WebGPUContext } from './context';
import { GPUShaderStage, GPUTextureUsage, GPUBufferUsage } from './compat';

export const depthPeelingShader = /* wgsl */ `
@group(3) @binding(0) var peelNear: texture_depth_2d;
@group(3) @binding(1) var peelFar: texture_depth_2d;
@group(3) @binding(2) var peelOpaque: texture_depth_2d;
@group(3) @binding(3) var peelFrontColor: texture_2d<f32>;
@group(3) @binding(4) var peelFrontEmission: texture_2d<f32>;
@group(3) @binding(5) var<uniform> peelSettings: vec4u;
fn peelSearch(depth: f32, pixel: vec2f) -> f32 {
    let p = vec2i(pixel);
    if (depth >= textureLoad(peelOpaque, p, 0)) { discard; }
    if (peelSettings.x > 0u && (depth <= textureLoad(peelNear, p, 0) || depth >= textureLoad(peelFar, p, 0))) { discard; }
    return depth;
}
struct PeelFragment { @location(0) color: vec4f, @location(1) emission: vec4f };
fn peelLayer(color: vec4f, emission: vec4f, depth: f32, pixel: vec2f, front: bool) -> PeelFragment {
    let p = vec2i(pixel); let nearest = textureLoad(peelNear, p, 0); let furthest = textureLoad(peelFar, p, 0);
    if (depth >= textureLoad(peelOpaque, p, 0)) { discard; }
    if (front) {
        if (depth != nearest) { discard; }
        let previous = textureLoad(peelFrontColor, p, 0); let previousEmission = textureLoad(peelFrontEmission, p, 0);
        return PeelFragment(previous + color * (1.0 - previous.a), previousEmission + emission * (1.0 - previousEmission.a));
    }
    if (depth != furthest || depth == nearest) { discard; }
    return PeelFragment(color, emission);
}
`;

const shader = /* wgsl */ `
@group(0) @binding(0) var frontColor: texture_2d<f32>;
@group(0) @binding(1) var frontEmission: texture_2d<f32>;
@group(0) @binding(2) var backColor: texture_2d<f32>;
@group(0) @binding(3) var backEmission: texture_2d<f32>;
@vertex fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
    let p = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    return vec4f(p[i], 0.0, 1.0);
}
struct Output { @location(0) color: vec4f, @location(1) emission: vec4f };
@fragment fn copy(@builtin(position) p: vec4f) -> Output {
    return Output(textureLoad(frontColor, vec2i(p.xy), 0), textureLoad(frontEmission, vec2i(p.xy), 0));
}
@fragment fn resolve(@builtin(position) p: vec4f) -> Output {
    let f = textureLoad(frontColor, vec2i(p.xy), 0); let e = textureLoad(frontEmission, vec2i(p.xy), 0);
    return Output(f + textureLoad(backColor, vec2i(p.xy), 0) * (1.0 - f.a),
        e + textureLoad(backEmission, vec2i(p.xy), 0) * (1.0 - e.a));
}
`;

/** Dual near/far peeling with full-precision depth and portable half-float MAX color blending. */
export class WebGPUDepthPeeling {
    readonly layout: GPUBindGroupLayout;
    readonly emptyLayout: GPUBindGroupLayout;
    readonly emptyGroup: GPUBindGroup;
    private readonly resolveLayout: GPUBindGroupLayout;
    private readonly copyPipeline: GPURenderPipeline;
    private readonly resolvePipelines: GPURenderPipeline[];
    private readonly settings: GPUBuffer[];
    private depths: GPUTexture[][] = [];
    private fronts: GPUTexture[][] = [];
    private backLayer: GPUTexture[] = [];
    private back: GPUTexture[] = [];
    static targets(): GPUColorTargetState[] {
        const blend: GPUBlendState = { color: { operation: 'max', srcFactor: 'one', dstFactor: 'one' }, alpha: { operation: 'max', srcFactor: 'one', dstFactor: 'one' } };
        return [0, 1].map(() => ({ format: 'rgba16float', blend }));
    }
    constructor(private readonly context: WebGPUContext) {
        const { device } = context;
        this.layout = device.createBindGroupLayout({ entries: [
            ...[0, 1, 2].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'depth' as const } })),
            ...[3, 4].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' as const } })),
            { binding: 5, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
        ] });
        this.emptyLayout = device.createBindGroupLayout({ entries: [] }); this.emptyGroup = device.createBindGroup({ layout: this.emptyLayout, entries: [] });
        this.resolveLayout = device.createBindGroupLayout({ entries: [0, 1, 2, 3].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' as const } })) });
        const module = device.createShaderModule({ label: 'molstar-depth-peeling-resolve', code: shader });
        const layout = device.createPipelineLayout({ bindGroupLayouts: [this.resolveLayout] });
        const blend: GPUBlendState = { color: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' }, alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' } };
        this.copyPipeline = device.createRenderPipeline({ layout, vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'copy', targets: [0, 1].map(() => ({ format: 'rgba16float', blend })) } });
        this.resolvePipelines = [15, 0, 8].map(mask => device.createRenderPipeline({ layout, vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'resolve', targets: [{ format: context.format, blend }, { format: 'rgba16float', blend, writeMask: mask }] } }));
        this.settings = [0, 1].map(initial => context.createBuffer(new Uint32Array([initial, 0, 0, 0]), GPUBufferUsage.UNIFORM));
    }
    private setSize(width: number, height: number) {
        if (this.depths[0]?.[0].width === width && this.depths[0][0].height === height) return;
        this.destroyTextures();
        const texture = (depth: boolean) => this.context.device.createTexture({ size: [width, height], format: depth ? 'depth32float' : 'rgba16float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC | GPUTextureUsage.COPY_DST });
        this.depths = [0, 1].map(() => [texture(true), texture(true)]);
        this.fronts = [0, 1].map(() => [texture(false), texture(false)]);
        this.backLayer = [texture(false), texture(false)]; this.back = [texture(false), texture(false)];
    }
    private viewport(pass: GPURenderPassEncoder, camera: Camera, height: number) {
        const v = camera.viewport; pass.setViewport(v.x, height - v.y - v.height, v.width, v.height, 0, 1);
    }
    private group(depths: GPUTexture[], opaque: GPUTexture, front: GPUTexture[], iteration: number) {
        return this.context.device.createBindGroup({ layout: this.layout, entries: [
            ...[...depths, opaque, ...front].map((texture, binding) => ({ binding, resource: texture.createView() })),
            { binding: 5, resource: { buffer: this.settings[iteration === 0 ? 0 : 1] } },
        ] });
    }
    private resolveGroup(front: GPUTexture[], back: GPUTexture[]) {
        return this.context.device.createBindGroup({ layout: this.resolveLayout, entries: [...front, ...back].map((texture, binding) => ({ binding, resource: texture.createView() })) });
    }
    render(encoder: GPUCommandEncoder, opaque: GPUTexture, color: GPUTexture, emission: GPUTexture, transparentColor: GPUTexture, camera: Camera, iterations: number, emitTransparent: boolean,
        draw: (pass: GPURenderPassEncoder, phase: 'near' | 'far' | 'front' | 'back', bindings: GPUBindGroup) => void) {
        this.setSize(color.width, color.height);
        const clear = (textures: GPUTexture[]) => { const pass = encoder.beginRenderPass({ colorAttachments: textures.map(texture => ({ view: texture.createView(), loadOp: 'clear' as const, storeOp: 'store' as const })) }); pass.end(); };
        clear(this.fronts[1]); clear(this.back);
        let front = this.fronts[1];
        for (let iteration = 0; iteration < iterations; iteration++) {
            const index = iteration % 2, previous = 1 - index;
            const search = this.group(this.depths[previous], opaque, front, iteration);
            for (const [slot, phase] of [[0, 'near'], [1, 'far']] as const) {
                const pass = encoder.beginRenderPass({ label: `molstar-peel-${phase}`, colorAttachments: [], depthStencilAttachment: {
                    view: this.depths[index][slot].createView(), depthClearValue: slot === 0 ? 1 : 0, depthLoadOp: 'clear', depthStoreOp: 'store',
                } });
                this.viewport(pass, camera, color.height); draw(pass, phase, search); pass.end();
            }
            const layer = this.group(this.depths[index], opaque, front, iteration);
            const next = this.fronts[index];
            for (let i = 0; i < 2; i++) encoder.copyTextureToTexture({ texture: front[i] }, { texture: next[i] }, [color.width, color.height]);
            for (const [textures, phase] of [[next, 'front'], [this.backLayer, 'back']] as const) {
                const pass = encoder.beginRenderPass({ label: `molstar-peel-${phase}`, colorAttachments: textures.map(texture => ({ view: texture.createView(), loadOp: phase === 'front' ? 'load' as const : 'clear' as const, storeOp: 'store' as const })) });
                this.viewport(pass, camera, color.height); draw(pass, phase, layer); pass.end();
            }
            const back = encoder.beginRenderPass({ label: 'molstar-peel-back-composite', colorAttachments: this.back.map(texture => ({ view: texture.createView(), loadOp: 'load', storeOp: 'store' })) });
            this.viewport(back, camera, color.height); back.setPipeline(this.copyPipeline); back.setBindGroup(0, this.resolveGroup(this.backLayer, this.backLayer)); back.draw(3); back.end();
            front = next;
        }
        const bindings = this.resolveGroup(front, this.back);
        for (const [output, transparentOnly] of [[transparentColor, true], [color, false]] as const) {
            const pass = encoder.beginRenderPass({ label: 'molstar-peel-resolve', colorAttachments: [output, emission].map((texture, i) => ({ view: texture.createView(), loadOp: transparentOnly && i === 0 ? 'clear' as const : 'load' as const, storeOp: 'store' as const })) });
            this.viewport(pass, camera, color.height); pass.setPipeline(this.resolvePipelines[transparentOnly ? 1 : emitTransparent ? 0 : 2]); pass.setBindGroup(0, bindings); pass.draw(3); pass.end();
        }
    }
    private destroyTextures() { for (const texture of [...this.depths.flat(), ...this.fronts.flat(), ...this.backLayer, ...this.back]) texture.destroy(); this.depths = []; this.fronts = []; this.backLayer = []; this.back = []; }
    dispose() { this.destroyTextures(); for (const buffer of this.settings) buffer.destroy(); }
}
