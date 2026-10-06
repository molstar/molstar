/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * SMAA 1x Medium port of Mol*'s three.js-derived SMAA v2.8 shaders.
 * MIT License Copyright (c) 2010-2020 three.js authors.
 */
import { Camera } from '../../mol-canvas3d/camera';
import { getSmaaLookupData } from './smaa-lookups';
import { SmaaProps } from '../../mol-canvas3d/passes/smaa';
import { GPUBufferUsage, GPUShaderStage, GPUTextureUsage } from './compat';
import { WebGPUContext } from './context';

export const smaaShader = /* wgsl */ `
struct Settings { dimensions: vec4f, viewport: vec4f };
@group(0) @binding(0) var<uniform> settings: Settings;
@group(0) @binding(1) var color: texture_2d<f32>;
@group(0) @binding(2) var edges: texture_2d<f32>;
@group(0) @binding(3) var weights: texture_2d<f32>;
@group(0) @binding(4) var area: texture_2d<f32>;
@group(0) @binding(5) var search: texture_2d<f32>;
@group(0) @binding(6) var linearSampler: sampler;
@group(0) @binding(7) var nearestSampler: sampler;
@vertex fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
    let positions = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    return vec4f(positions[i], 0.0, 1.0);
}
fn inside(p: vec2f) -> bool { return all(p >= settings.viewport.xy) && all(p < settings.viewport.xy + settings.viewport.zw); }
// Keep SMAA's original bottom-left search coordinates; texture attachments are top-left.
fn textureUv(uv: vec2f) -> vec2f {
    let p = vec2f(uv.x, 1.0 - uv.y);
    return clamp(p, (settings.viewport.xy + 0.5) / settings.dimensions.xy, (settings.viewport.xy + settings.viewport.zw - 0.5) / settings.dimensions.xy);
}
fn sampleColor(uv: vec2f) -> vec4f { return textureSampleLevel(color, linearSampler, textureUv(uv), 0.0); }
fn sampleEdges(uv: vec2f) -> vec2f { return textureSampleLevel(edges, linearSampler, textureUv(uv), 0.0).rg; }
fn sampleWeights(uv: vec2f) -> vec4f { return textureSampleLevel(weights, linearSampler, textureUv(uv), 0.0); }
fn maxChannel(v: vec3f) -> f32 { return max(v.r, max(v.g, v.b)); }
fn texcoord(p: vec2f) -> vec2f { return vec2f(p.x / settings.dimensions.x, 1.0 - p.y / settings.dimensions.y); }
@fragment fn detectEdges(@builtin(position) p: vec4f) -> @location(0) vec4f {
    if (!inside(p.xy)) { return vec4f(0.0); }
    let uv = texcoord(p.xy); let pixel = 1.0 / settings.dimensions.xy;
    let c = sampleColor(uv).rgb;
    let delta = vec4f(maxChannel(abs(c - sampleColor(uv - vec2f(pixel.x, 0.0)).rgb)),
        maxChannel(abs(c - sampleColor(uv + vec2f(0.0, pixel.y)).rgb)),
        maxChannel(abs(c - sampleColor(uv + vec2f(pixel.x, 0.0)).rgb)),
        maxChannel(abs(c - sampleColor(uv - vec2f(0.0, pixel.y)).rgb)));
    let initial = step(vec2f(settings.dimensions.z), delta.xy);
    if (dot(initial, vec2f(1.0)) == 0.0) { return vec4f(0.0); }
    let farDelta = vec2f(maxChannel(abs(c - sampleColor(uv - vec2f(2.0 * pixel.x, 0.0)).rgb)), maxChannel(abs(c - sampleColor(uv + vec2f(0.0, 2.0 * pixel.y)).rgb)));
    let maximum = max(max(delta.x, delta.y), max(max(delta.z, delta.w), max(farDelta.x, farDelta.y)));
    return vec4f(initial * step(vec2f(0.5 * maximum), delta.xy), 0.0, 0.0);
}
fn searchLength(e: vec2f, bias: f32) -> f32 {
    return 255.0 * textureSampleLevel(search, nearestSampler, vec2f(bias + e.x * 0.5, e.y), 0.0).r;
}
fn searchX(start: vec2f, end: f32, right: bool) -> f32 {
    let pixel = 1.0 / settings.dimensions.xy;
    let sign = select(-1.0, 1.0, right);
    var uv = start; var e = vec2f(0.0, 1.0);
    for (var i = 0u; i < u32(settings.dimensions.w); i++) {
        e = sampleEdges(uv); uv.x += sign * 2.0 * pixel.x;
        let searching = select(uv.x > end, uv.x < end, right);
        if (!(searching && e.y > 0.8281 && e.x == 0.0)) { break; }
    }
    return uv.x - sign * (3.25 - searchLength(e, select(0.0, 0.5, right))) * pixel.x;
}
fn searchY(start: vec2f, end: f32, down: bool) -> f32 {
    let pixel = 1.0 / settings.dimensions.xy;
    let sign = select(1.0, -1.0, down);
    var uv = start; var e = vec2f(1.0, 0.0);
    for (var i = 0u; i < u32(settings.dimensions.w); i++) {
        e = sampleEdges(uv); uv.y += sign * 2.0 * pixel.y;
        let searching = select(uv.y > end, uv.y < end, down);
        if (!(searching && e.x > 0.8281 && e.y == 0.0)) { break; }
    }
    return uv.y - sign * (3.25 - searchLength(e.yx, select(0.0, 0.5, down))) * pixel.y;
}
fn areaWeights(dist: vec2f, e1: f32, e2: f32) -> vec2f {
    let uv = (16.0 * round(4.0 * vec2f(e1, e2)) + dist + 0.5) / vec2f(160.0, 560.0);
    return textureSampleLevel(area, linearSampler, uv, 0.0).rg;
}
@fragment fn calculateWeights(@builtin(position) p: vec4f) -> @location(0) vec4f {
    if (!inside(p.xy)) { return vec4f(0.0); }
    let uv = texcoord(p.xy); let pixel = 1.0 / settings.dimensions.xy;
    let o0 = uv.xyxy + pixel.xyxy * vec4f(-0.25, 0.125, 1.25, 0.125);
    let o1 = uv.xyxy + pixel.xyxy * vec4f(-0.125, 0.25, -0.125, -1.25);
    let ends = vec4f(o0.xz, o1.yw) + vec4f(-2.0, 2.0, -2.0, 2.0) * pixel.xxyy * settings.dimensions.w;
    let e = sampleEdges(uv); var result = vec4f(0.0);
    if (e.y > 0.0) {
        let left = searchX(o0.xy, ends.x, false); let right = searchX(o0.zw, ends.y, true);
        let e1 = sampleEdges(vec2f(left, o1.y)).x;
        let e2 = sampleEdges(vec2f(right + pixel.x, o1.y - pixel.y)).x;
        let dist = sqrt(abs(vec2f(left, right) / pixel.x - uv.x / pixel.x));
        result = vec4f(areaWeights(dist, e1, e2), result.ba);
    }
    if (e.x > 0.0) {
        let up = searchY(o1.xy, ends.z, false); let down = searchY(o1.zw, ends.w, true);
        let e1 = sampleEdges(vec2f(o0.x, up)).y;
        let e2 = sampleEdges(vec2f(o0.x, down)).y;
        let dist = sqrt(abs(vec2f(up, down) / pixel.y - uv.y / pixel.y));
        result = vec4f(result.rg, areaWeights(dist, e1, e2));
    }
    return result;
}
@fragment fn blend(@builtin(position) p: vec4f) -> @location(0) vec4f {
    let original = textureLoad(color, vec2i(p.xy), 0);
    if (!inside(p.xy)) { return original; }
    let uv = texcoord(p.xy); let pixel = 1.0 / settings.dimensions.xy;
    let currentWeights = sampleWeights(uv);
    let a = vec4f(currentWeights.x, sampleWeights(uv - vec2f(0.0, pixel.y)).g, currentWeights.z, sampleWeights(uv + vec2f(pixel.x, 0.0)).a);
    if (dot(a, vec4f(1.0)) < 0.00001) { return original; }
    var offset = vec2f(select(-a.b, a.a, a.a > a.b), select(a.r, -a.g, a.g > a.r));
    if (abs(offset.x) > abs(offset.y)) { offset.y = 0.0; } else { offset.x = 0.0; }
    let opposite = sampleColor(uv + sign(offset) * pixel);
    let weight = max(abs(offset.x), abs(offset.y));
    // Gamma-correct premultiplied blending also preserves transparent silhouette coverage.
    let c = pow(original.rgb / max(original.a, 0.000001), vec3f(2.2)) * original.a;
    let o = pow(opposite.rgb / max(opposite.a, 0.000001), vec3f(2.2)) * opposite.a;
    let alpha = mix(original.a, opposite.a, weight);
    let rgb = pow(max(mix(c, o, weight) / max(alpha, 0.000001), vec3f(0.0)), vec3f(1.0 / 2.2)) * alpha;
    return vec4f(rgb, alpha);
}
`;

interface LookupTextures { area: GPUTexture, search: GPUTexture }
interface LookupEntry { textures: Promise<LookupTextures>, references: number }
const lookupCache = new WeakMap<WebGPUContext, LookupEntry>();

async function loadLookups(context: WebGPUContext): Promise<LookupTextures> {
    const data = await getSmaaLookupData();
    const load = (image: { width: number, height: number, array: Uint8Array }, label: string) => {
        const texture = context.device.createTexture({ label, size: { width: image.width, height: image.height }, format: 'rgba8unorm', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST });
        try {
            context.device.queue.writeTexture({ texture }, image.array, { bytesPerRow: image.width * 4, rowsPerImage: image.height }, { width: image.width, height: image.height });
            return texture;
        } catch (error) { texture.destroy(); throw error; }
    };
    const area = load(data.area, 'molstar-smaa-area');
    try { return { area, search: load(data.search, 'molstar-smaa-search') }; } catch (error) { area.destroy(); throw error; }
}

/** Three native SMAA render passes with shared lookup tables, no WebGL resources. */
export class WebGPUSmaa {
    private readonly layout: GPUBindGroupLayout;
    private readonly pipelines: Record<'detectEdges' | 'calculateWeights' | 'blend', GPURenderPipeline>;
    private readonly linear: GPUSampler;
    private readonly nearest: GPUSampler;
    private readonly uniforms: GPUBuffer;
    private edges?: GPUTexture;
    private weights?: GPUTexture;
    private output?: GPUTexture;
    private disposed = false;

    private constructor(private readonly context: WebGPUContext, shader: GPUShaderModule, private readonly area: GPUTexture, private readonly search: GPUTexture, private readonly releaseLookups: () => void) {
        const { device } = context;
        this.layout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            ...[1, 2, 3, 4, 5].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' as const } })),
            ...[6, 7].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' as const } })),
        ] });
        const layout = device.createPipelineLayout({ bindGroupLayouts: [this.layout] });
        const pipeline = (entryPoint: string, format: GPUTextureFormat) => device.createRenderPipeline({ label: `molstar-smaa-${entryPoint}`, layout, vertex: { module: shader, entryPoint: 'vs' }, fragment: { module: shader, entryPoint, targets: [{ format }] }, primitive: { topology: 'triangle-list' } });
        this.pipelines = { detectEdges: pipeline('detectEdges', 'rgba8unorm'), calculateWeights: pipeline('calculateWeights', 'rgba8unorm'), blend: pipeline('blend', context.format) };
        this.uniforms = device.createBuffer({ size: 32, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.linear = device.createSampler({ minFilter: 'linear', magFilter: 'linear' });
        this.nearest = device.createSampler({ minFilter: 'nearest', magFilter: 'nearest' });
    }

    static async create(context: WebGPUContext, shader: GPUShaderModule) {
        let entry = lookupCache.get(context);
        if (!entry) { entry = { textures: loadLookups(context), references: 0 }; lookupCache.set(context, entry); }
        const shared = entry;
        shared.references++;
        const release = () => {
            if (--shared.references === 0) {
                lookupCache.delete(context);
                shared.textures.then(({ area, search }) => { area.destroy(); search.destroy(); }).catch(() => {});
            }
        };
        try {
            const { area, search } = await shared.textures;
            return new WebGPUSmaa(context, shader, area, search, release);
        } catch (error) { release(); throw error; }
    }

    render(encoder: GPUCommandEncoder, color: GPUTexture, camera: Camera, props: SmaaProps) {
        if (this.disposed) throw new Error('WebGPU SMAA has been disposed.');
        if (this.output?.width !== color.width || this.output.height !== color.height) {
            this.destroyTargets();
            const target = (format: GPUTextureFormat) => this.context.device.createTexture({ size: { width: color.width, height: color.height }, format, usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC });
            this.edges = target('rgba8unorm'); this.weights = target('rgba8unorm'); this.output = target(this.context.format);
        }
        const v = camera.viewport;
        this.context.device.queue.writeBuffer(this.uniforms, 0, new Float32Array([color.width, color.height, props.edgeThreshold, props.maxSearchSteps, v.x, color.height - v.y - v.height, v.width, v.height]));
        const run = (entry: 'detectEdges' | 'calculateWeights' | 'blend', target: GPUTexture) => {
            const bindings = this.context.device.createBindGroup({ layout: this.layout, entries: [
                { binding: 0, resource: { buffer: this.uniforms } }, { binding: 1, resource: color.createView() },
                { binding: 2, resource: (entry === 'calculateWeights' ? this.edges! : color).createView() },
                { binding: 3, resource: (entry === 'blend' ? this.weights! : color).createView() },
                { binding: 4, resource: this.area.createView() }, { binding: 5, resource: this.search.createView() },
                { binding: 6, resource: this.linear }, { binding: 7, resource: this.nearest },
            ] });
            const pass = encoder.beginRenderPass({ label: `molstar-smaa-${entry}`, colorAttachments: [{ view: target.createView(), loadOp: 'clear', storeOp: 'store' }] });
            pass.setPipeline(this.pipelines[entry]); pass.setBindGroup(0, bindings); pass.draw(3); pass.end();
        };
        run('detectEdges', this.edges!); run('calculateWeights', this.weights!); run('blend', this.output!);
        return this.output!;
    }

    private destroyTargets() { this.edges?.destroy(); this.weights?.destroy(); this.output?.destroy(); }
    dispose() { if (this.disposed) return; this.disposed = true; this.destroyTargets(); this.uniforms.destroy(); this.releaseLookups(); }
}
