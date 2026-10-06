/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { getCanvasModule } from '../../mol-geo/geometry/text/font-atlas';
import { Subject } from 'rxjs';
import { BackgroundProps } from '../../mol-canvas3d/passes/background';
import { Camera } from '../../mol-canvas3d/camera';
import { Asset, AssetManager } from '../../mol-util/assets';
import { Color } from '../../mol-util/color';
import { Mat3, Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Euler } from '../../mol-math/linear-algebra/3d/euler';
import { degToRad } from '../../mol-math/misc';
import { WebGPUContext } from './context';
import { GPUBufferUsage, GPUShaderStage, GPUTextureUsage } from './compat';

const shader = /* wgsl */ `
struct Settings { dimensions: vec4f, viewport: vec4f, first: vec4f, second: vec4f, options: vec4f, background: vec4f, inverse: mat4x4f, rotation: mat4x4f };
@group(0) @binding(0) var<uniform> settings: Settings;
@group(0) @binding(1) var scene: texture_2d<f32>;
@group(0) @binding(2) var image: texture_2d<f32>;
@group(0) @binding(3) var sky: texture_cube<f32>;
@group(0) @binding(4) var linearSampler: sampler;
@vertex fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
    let positions = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    return vec4f(positions[i], 0.0, 1.0);
}
@fragment fn fs(@builtin(position) position: vec4f) -> @location(0) vec4f {
    let original = textureLoad(scene, vec2i(position.xy), 0);
    let solid = vec4f(settings.background.rgb * settings.background.a, settings.background.a);
    if (any(position.xy < settings.viewport.xy) || any(position.xy >= settings.viewport.xy + settings.viewport.zw)) { return original + solid * (1.0 - original.a); }
    let uv = select(position.xy / settings.dimensions.xy, (position.xy - settings.viewport.xy) / settings.viewport.zw, settings.dimensions.z > 0.0);
    var rgb = vec3f(0.0); var alpha = settings.options.z;
    if (settings.options.x < 3.0) {
        var d = 1.0 - uv.y + 1.0 - settings.first.w * 2.0;
        if (settings.options.x == 2.0) { d = 1.0 - (distance(vec2f(0.5), uv) + settings.first.w - 0.5); }
        rgb = mix(settings.second.rgb, settings.first.rgb, clamp(d, 0.0, 1.0));
        let value = dot(vec2f(171.0, 231.0), vec2f(position.x, settings.dimensions.y - position.y));
        rgb += fract(vec3f(value) / vec3f(103.0, 71.0, 97.0)) / 255.0;
        alpha = 1.0;
    } else {
        if (settings.options.x == 3.0) {
            let aspect = settings.dimensions.x / settings.dimensions.y;
            let imageAspect = f32(textureDimensions(image).x) / f32(textureDimensions(image).y);
            let scale = select(vec2f(1.0, imageAspect / aspect), vec2f(aspect / imageAspect, 1.0), aspect < imageAspect);
            rgb = textureSampleLevel(image, linearSampler, uv * scale + (1.0 - scale) * 0.5, settings.options.y * 8.0).rgb;
        } else {
            let local = (position.xy - settings.viewport.xy) / settings.viewport.zw;
            let point = settings.inverse * vec4f(local.x * 2.0 - 1.0, 1.0 - local.y * 2.0, 1.0, 1.0);
            let direction = (settings.rotation * vec4f(normalize(point.xyz / point.w), 0.0)).xyz;
            rgb = textureSampleLevel(sky, linearSampler, direction, settings.options.y * 8.0).rgb;
        }
        let intensity = vec3f(dot(rgb, vec3f(0.2125, 0.7154, 0.0721)));
        rgb = mix(intensity, rgb, 1.0 + settings.options.w) + settings.second.w;
    }
    let environment = vec4f(clamp(rgb, vec3f(0.0), vec3f(1.0)) * alpha, alpha) + solid * (1.0 - alpha);
    return original + environment * (1.0 - original.a);
}
`;
const mipShader = /* wgsl */ `
@group(0) @binding(0) var source: texture_2d<f32>;
@group(0) @binding(1) var linearSampler: sampler;
struct Vertex { @builtin(position) position: vec4f, @location(0) uv: vec2f };
@vertex fn vs(@builtin(vertex_index) i: u32) -> Vertex {
    let vertices = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    var out: Vertex; out.position = vec4f(vertices[i], 0.0, 1.0); out.uv = vec2f(vertices[i].x * 0.5 + 0.5, 0.5 - vertices[i].y * 0.5); return out;
}
@fragment fn fs(in: Vertex) -> @location(0) vec4f { return textureSampleLevel(source, linearSampler, in.uv, 0.0); }
`;

async function decodeImage(data: Uint8Array<ArrayBuffer>): Promise<{ width: number, height: number, bitmap?: ImageBitmap, pixels?: Uint8ClampedArray, close(): void }> {
    if (typeof createImageBitmap !== 'undefined') {
        const bitmap = await createImageBitmap(new Blob([data]), { premultiplyAlpha: 'none', colorSpaceConversion: 'none' });
        return { width: bitmap.width, height: bitmap.height, bitmap, close: () => bitmap.close() };
    }
    const module = getCanvasModule();
    if (!module.loadImage) throw new Error('Headless background assets require a canvas module with loadImage.');
    const image = await module.loadImage(data), canvas = module.createCanvas(image.width, image.height);
    const context = canvas.getContext('2d');
    context.drawImage(image, 0, 0);
    const pixels: Uint8ClampedArray = context.getImageData(0, 0, image.width, image.height).data;
    return { width: image.width, height: image.height, pixels, close: () => { } };
}

/** Native gradient, image and cubemap backgrounds; assets never acquire WebGL. */
export class WebGPUBackground {
    readonly changed = new Subject<void>();
    private readonly assets: AssetManager;
    private readonly ownsAssets: boolean;
    private readonly layout: GPUBindGroupLayout;
    private readonly pipeline: GPURenderPipeline;
    private readonly mipPipeline: GPURenderPipeline;
    private readonly sampler: GPUSampler;
    private readonly uniform: GPUBuffer;
    private readonly emptyImage: GPUTexture;
    private readonly emptySky: GPUTexture;
    private color?: GPUTexture;
    private texture?: GPUTexture;
    private wrappers: Asset.Wrapper<'binary'>[] = [];
    private key = '';
    private revision = 0;
    private pending?: Promise<void>;
    private disposed = false;
    constructor(private readonly context: WebGPUContext, assets?: AssetManager) {
        this.assets = assets ?? new AssetManager(); this.ownsAssets = !assets;
        const { device } = context;
        this.layout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            ...[1, 2].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' as const } })),
            { binding: 3, visibility: GPUShaderStage.FRAGMENT, texture: { viewDimension: 'cube', sampleType: 'float' } },
            { binding: 4, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
        ] });
        const module = device.createShaderModule({ label: 'molstar-background', code: shader });
        this.pipeline = device.createRenderPipeline({ layout: device.createPipelineLayout({ bindGroupLayouts: [this.layout] }), vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'fs', targets: [{ format: context.format }] } });
        const mipModule = device.createShaderModule({ label: 'molstar-background-mips', code: mipShader });
        this.mipPipeline = device.createRenderPipeline({ layout: 'auto', vertex: { module: mipModule, entryPoint: 'vs' }, fragment: { module: mipModule, entryPoint: 'fs', targets: [{ format: 'rgba8unorm' }] } });
        this.sampler = device.createSampler({ minFilter: 'linear', magFilter: 'linear', mipmapFilter: 'linear' });
        this.uniform = device.createBuffer({ size: 224, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.emptyImage = device.createTexture({ size: [1, 1], format: 'rgba8unorm', usage: GPUTextureUsage.TEXTURE_BINDING });
        this.emptySky = device.createTexture({ size: [1, 1, 6], format: 'rgba8unorm', usage: GPUTextureUsage.TEXTURE_BINDING });
    }
    private clear() {
        this.revision++; this.texture?.destroy(); this.texture = undefined;
        for (const wrapper of this.wrappers) wrapper.dispose(); this.wrappers = [];
        this.key = ''; this.pending = undefined;
    }
    private async load(assets: Asset[], cube: boolean, revision: number) {
        const settled = await Promise.allSettled(assets.map(async asset => {
            const wrapper = await this.assets.resolve(asset, 'binary').run();
            try { return { wrapper, image: await decodeImage(wrapper.data) }; } catch (error) { wrapper.dispose(); throw error; }
        }));
        const loaded = settled.flatMap(result => result.status === 'fulfilled' ? [result.value] : []);
        let texture: GPUTexture | undefined;
        try {
            const failed = settled.find(result => result.status === 'rejected');
            if (failed?.status === 'rejected') throw failed.reason;
            if (this.disposed || revision !== this.revision) return;
            const { width, height } = loaded[0].image;
            if (width > this.context.device.limits.maxTextureDimension2D || height > this.context.device.limits.maxTextureDimension2D) throw new Error('WebGPU background exceeds the device texture dimension limit.');
            if (cube && (width !== height || loaded.some(item => item.image.width !== width || item.image.height !== height))) throw new Error('WebGPU skybox faces must be equally sized square images.');
            const levels = Math.floor(Math.log2(Math.max(width, height))) + 1, { device } = this.context;
            texture = device.createTexture({ label: 'molstar-background-asset', size: [width, height, loaded.length], mipLevelCount: levels, format: 'rgba8unorm', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.RENDER_ATTACHMENT });
            for (let layer = 0; layer < loaded.length; layer++) {
                const image = loaded[layer].image;
                if (image.bitmap) device.queue.copyExternalImageToTexture({ source: image.bitmap }, { texture, origin: [0, 0, layer], premultipliedAlpha: false }, [width, height]);
                else device.queue.writeTexture({ texture, origin: [0, 0, layer] }, image.pixels!, { bytesPerRow: width * 4 }, [width, height]);
            }
            const encoder = device.createCommandEncoder({ label: 'molstar-background-mips' });
            for (let layer = 0; layer < loaded.length; layer++) for (let level = 1; level < levels; level++) {
                const bindings = device.createBindGroup({ layout: this.mipPipeline.getBindGroupLayout(0), entries: [
                    { binding: 0, resource: texture.createView({ dimension: '2d', baseArrayLayer: layer, arrayLayerCount: 1, baseMipLevel: level - 1, mipLevelCount: 1 }) }, { binding: 1, resource: this.sampler },
                ] });
                const pass = encoder.beginRenderPass({ colorAttachments: [{ view: texture.createView({ dimension: '2d', baseArrayLayer: layer, arrayLayerCount: 1, baseMipLevel: level, mipLevelCount: 1 }), loadOp: 'clear', storeOp: 'store' }] });
                pass.setPipeline(this.mipPipeline); pass.setBindGroup(0, bindings); pass.draw(3); pass.end();
            }
            device.queue.submit([encoder.finish()]);
            this.texture = texture; texture = undefined;
            this.wrappers = loaded.map(item => item.wrapper);
            this.changed.next();
        } finally {
            texture?.destroy();
            for (const item of loaded) { item.image.close(); if (!this.wrappers.includes(item.wrapper)) item.wrapper.dispose(); }
        }
    }
    update(props: BackgroundProps) {
        const variant = props.variant;
        if (variant.name === 'off' || variant.name === 'horizontalGradient' || variant.name === 'radialGradient') {
            if (this.key) this.clear();
            return variant.name !== 'off';
        }
        let assets: Asset[] = [], key = '';
        if (variant.name === 'image') {
            const source = variant.params.source;
            if (source.params) {
                const asset = source.name === 'url' ? this.assets.tryFindFilename(source.params) ?? Asset.getUrlAsset(this.assets, source.params) : source.params;
                assets = [asset]; key = `image/${source.name}/${source.name === 'url' ? source.params : asset.id}`;
            }
        } else {
            const faces = variant.params.faces;
            // WebGPU cube layer order: https://www.w3.org/TR/webgpu/#dom-gputextureviewdimension-cube
            const ordered = [faces.params.px, faces.params.nx, faces.params.py, faces.params.ny, faces.params.pz, faces.params.nz];
            if (ordered.every(Boolean)) {
                assets = faces.name === 'urls' ? (ordered as string[]).map(url => Asset.getUrlAsset(this.assets, url)) : ordered as Asset[];
                key = `sky/${faces.name}/${ordered.map(face => typeof face === 'string' ? face : face!.id).join('|')}`;
            }
        }
        if (!key) { if (this.key) this.clear(); return false; }
        if (this.key !== key) {
            this.clear(); this.key = key;
            this.pending = this.load(assets, variant.name === 'skybox', this.revision);
            this.pending.catch(() => { if (!this.disposed) this.changed.next(); });
        }
        return !!this.texture;
    }
    async ready(props: BackgroundProps) { this.update(props); await this.pending; }
    render(encoder: GPUCommandEncoder, source: GPUTexture, camera: Camera, props: BackgroundProps, transparent: boolean, background: Color) {
        if (!this.update(props)) return source;
        if (!this.color || this.color.width !== source.width || this.color.height !== source.height) {
            this.color?.destroy(); this.color = this.context.device.createTexture({ size: [source.width, source.height], format: this.context.format, usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC });
        }
        const data = new Float32Array(56), variant = props.variant, v = camera.viewport;
        data.set([source.width, source.height, variant.name !== 'skybox' && variant.name !== 'off' && variant.params.coverage === 'viewport' ? 1 : 0, 0]);
        data.set([v.x, source.height - v.y - v.height, v.width, v.height], 4);
        data.set([...Color.toRgbNormalized(background), transparent ? 0 : 1], 20);
        if (variant.name === 'horizontalGradient' || variant.name === 'radialGradient') {
            const p = variant.params;
            data[16] = variant.name === 'horizontalGradient' ? 1 : 2;
            data.set([...Color.toRgbNormalized('topColor' in p ? p.topColor : p.centerColor), p.ratio], 8);
            data.set(Color.toRgbNormalized('bottomColor' in p ? p.bottomColor : p.edgeColor), 12);
        } else if (variant.name === 'image' || variant.name === 'skybox') {
            const p = variant.params;
            data.set([variant.name === 'image' ? 3 : 4, p.blur, p.opacity, p.saturation], 16); data[15] = p.lightness;
            if (variant.name === 'skybox') {
                const perspective = camera.state.mode === 'orthographic' ? new Camera({ ...camera.state, mode: 'perspective' }, camera.viewport) : camera;
                perspective.update();
                const view = Mat4();
                if (Mat4.isZero(camera.headRotation)) {
                    const direction = Vec3.setMagnitude(Vec3(), Vec3.sub(Vec3(), perspective.state.position, perspective.state.target), 0.1);
                    Mat4.lookAt(view, direction, Vec3(), perspective.state.up);
                } else Mat4.invert(view, camera.headRotation);
                data.set(Mat4.invert(Mat4(), Mat4.mul(Mat4(), perspective.projection, view)), 24);
                const r = variant.params.rotation;
                const rotation = Mat3.fromEuler(Mat3(), Euler.create(degToRad(r.x), degToRad(r.y), degToRad(r.z)), 'XYZ');
                const matrix = Mat4.identity(); for (let column = 0; column < 3; column++) for (let row = 0; row < 3; row++) matrix[column * 4 + row] = rotation[column * 3 + row];
                data.set(matrix, 40);
            }
        }
        this.context.device.queue.writeBuffer(this.uniform, 0, data);
        const bindings = this.context.device.createBindGroup({ layout: this.layout, entries: [
            { binding: 0, resource: { buffer: this.uniform } }, { binding: 1, resource: source.createView() },
            { binding: 2, resource: (variant.name === 'image' ? this.texture! : this.emptyImage).createView() },
            { binding: 3, resource: (variant.name === 'skybox' ? this.texture! : this.emptySky).createView({ dimension: 'cube' }) }, { binding: 4, resource: this.sampler },
        ] });
        const pass = encoder.beginRenderPass({ label: 'molstar-background-compose', colorAttachments: [{ view: this.color.createView(), loadOp: 'clear', storeOp: 'store' }] });
        pass.setPipeline(this.pipeline); pass.setBindGroup(0, bindings); pass.draw(3); pass.end();
        return this.color;
    }
    dispose() { this.disposed = true; this.clear(); this.color?.destroy(); this.emptyImage.destroy(); this.emptySky.destroy(); this.uniform.destroy(); this.changed.complete(); if (this.ownsAssets) this.assets.dispose(); }
}
