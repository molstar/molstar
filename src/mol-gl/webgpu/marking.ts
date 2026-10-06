/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { MarkingProps } from '../../mol-canvas3d/passes/marking';
import { Camera } from '../../mol-canvas3d/camera';
import { Color } from '../../mol-util/color';
import { WebGPUContext } from './context';
import { GPUBufferUsage, GPUShaderStage, GPUTextureUsage } from './compat';

const shader = /* wgsl */ `
struct Settings { dimensions: vec4f, viewport: vec4f, highlight: vec4f, selection: vec4f, inner: vec4f };
@group(0) @binding(0) var<uniform> settings: Settings;
@group(0) @binding(1) var source: texture_2d<f32>;
@group(0) @binding(2) var mask: texture_2d<f32>;
@vertex fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
    let positions = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    return vec4f(positions[i], 0.0, 1.0);
}
fn inside(p: vec2f) -> bool { return all(p >= settings.viewport.xy) && all(p < settings.viewport.xy + settings.viewport.zw); }
fn sampleMask(p: vec2i) -> vec4f { return textureLoad(mask, clamp(p, vec2i(settings.viewport.xy), vec2i(settings.viewport.xy + settings.viewport.zw) - 1), 0); }
@fragment fn edge(@builtin(position) position: vec4f) -> @location(0) vec4f {
    if (!inside(position.xy)) { return vec4f(0.0); }
    let p = vec2i(position.xy); let scale = i32(settings.dimensions.z);
    let c0 = sampleMask(p); let c1 = sampleMask(p + vec2i(scale, 0)); let c2 = sampleMask(p - vec2i(scale, 0));
    let c3 = sampleMask(p + vec2i(0, scale)); let c4 = sampleMask(p - vec2i(0, scale));
    if (length(vec2f(c1.r - c2.r, c3.r - c4.r)) == 0.0) { return vec4f(0.0); }
    return vec4f(select(0.0, 1.0, min(min(c1.g, c2.g), min(c3.g, c4.g)) > 0.001), c0.r,
        min(min(c1.b, c2.b), min(c3.b, c4.b)), min(min(c1.a, c2.a), min(c3.a, c4.a)));
}
@fragment fn overlay(@builtin(position) position: vec4f) -> @location(0) vec4f {
    let p = vec2i(position.xy); let original = textureLoad(source, p, 0);
    if (!inside(position.xy)) { return original; }
    let edge = textureLoad(mask, p, 0);
    let tint = select(settings.selection, settings.highlight, edge.b == 1.0);
    let rgb = tint.rgb * select(settings.inner.x, 1.0, edge.g > 0.0);
    let alpha = clamp(select(1.0, settings.dimensions.w, edge.r == 1.0) * edge.a * tint.a, 0.0, 1.0);
    return vec4f(rgb * alpha, alpha) + original * (1.0 - alpha);
}
`;

/** Primitive selection/hover masks and fog-aware visible/hidden marking edges. */
export class WebGPUMarking {
    readonly geometryLayout: GPUBindGroupLayout;
    private readonly layout: GPUBindGroupLayout;
    private readonly uniform: GPUBuffer;
    private readonly edge: GPURenderPipeline;
    private readonly overlay: GPURenderPipeline;
    depth?: GPUTexture;
    maskDepth?: GPUTexture;
    mask?: GPUTexture;
    private edges?: GPUTexture;
    private color?: GPUTexture;
    constructor(private readonly context: WebGPUContext) {
        const { device } = context;
        this.geometryLayout = device.createBindGroupLayout({ entries: [{ binding: 0, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'depth' } }] });
        this.layout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            ...[1, 2].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' as const } })),
        ] });
        this.uniform = device.createBuffer({ size: 80, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        const module = device.createShaderModule({ label: 'molstar-marking', code: shader });
        const layout = device.createPipelineLayout({ bindGroupLayouts: [this.layout] });
        const pipeline = (entryPoint: string, format: GPUTextureFormat) => device.createRenderPipeline({ layout, vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint, targets: [{ format }] } });
        this.edge = pipeline('edge', 'rgba16float'); this.overlay = pipeline('overlay', context.format);
    }
    setSize(width: number, height: number) {
        if (this.depth?.width === width && this.depth.height === height) return;
        this.destroyTextures();
        const size = { width, height }, usage = GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING;
        this.depth = this.context.device.createTexture({ label: 'molstar-unmarked-depth', size, format: 'depth32float', usage });
        this.maskDepth = this.context.device.createTexture({ label: 'molstar-marked-depth', size, format: 'depth32float', usage: GPUTextureUsage.RENDER_ATTACHMENT });
        this.mask = this.context.device.createTexture({ label: 'molstar-marking-mask', size, format: 'rgba16float', usage });
        this.edges = this.context.device.createTexture({ label: 'molstar-marking-edges', size, format: 'rgba16float', usage });
        this.color = this.context.device.createTexture({ label: 'molstar-marking-color', size, format: this.context.format, usage: usage | GPUTextureUsage.COPY_SRC });
    }
    geometryBindings() { return this.context.device.createBindGroup({ layout: this.geometryLayout, entries: [{ binding: 0, resource: this.depth!.createView() }] }); }
    render(encoder: GPUCommandEncoder, source: GPUTexture, camera: Camera, props: MarkingProps, pixelRatio: number) {
        const settings = new Float32Array(20), v = camera.viewport;
        settings.set([source.width, source.height, Math.max(1, Math.round(props.edgeScale * pixelRatio)), props.ghostEdgeStrength]);
        settings.set([v.x, source.height - v.y - v.height, v.width, v.height], 4);
        settings.set([...Color.toRgbNormalized(props.highlightEdgeColor), props.highlightEdgeStrength], 8);
        settings.set([...Color.toRgbNormalized(props.selectEdgeColor), props.selectEdgeStrength], 12);
        settings[16] = props.innerEdgeFactor;
        this.context.device.queue.writeBuffer(this.uniform, 0, settings);
        const run = (name: string, pipeline: GPURenderPipeline, mask: GPUTexture, output: GPUTexture) => {
            const bindings = this.context.device.createBindGroup({ layout: this.layout, entries: [
                { binding: 0, resource: { buffer: this.uniform } }, { binding: 1, resource: source.createView() }, { binding: 2, resource: mask.createView() },
            ] });
            const pass = encoder.beginRenderPass({ label: `molstar-marking-${name}`, colorAttachments: [{ view: output.createView(), loadOp: 'clear', storeOp: 'store' }] });
            pass.setPipeline(pipeline); pass.setBindGroup(0, bindings); pass.draw(3); pass.end();
        };
        run('edges', this.edge, this.mask!, this.edges!); run('overlay', this.overlay, this.edges!, this.color!);
        return this.color!;
    }
    private destroyTextures() { for (const t of [this.depth, this.maskDepth, this.mask, this.edges, this.color]) t?.destroy(); }
    dispose() { this.destroyTextures(); this.uniform.destroy(); }
}
