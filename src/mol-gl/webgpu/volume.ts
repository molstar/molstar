/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { getWebGPUModelView } from './camera';
import { Camera } from '../../mol-canvas3d/camera';
import { Mat4 } from '../../mol-math/linear-algebra';
import { fromHalfFloat } from '../../mol-util/number-conversion';
import { GraphicsRenderObject } from '../render-object';
import { RenderableValues } from '../renderable/schema';
import { GPUBufferUsage, GPUShaderStage, GPUTextureUsage } from './compat';
import { WebGPUContext } from './context';
import { value } from './geometry';
import { createClipData } from './clip';
import { WebGPUTextureData } from './texture-data';
import { WebGPUWeightedTransparency } from './transparency';
import { WebGPUDepthPeeling } from './depth-peeling';

interface VolumeItem {
    object: GraphicsRenderObject
    buffers: GPUBuffer[]
    textures: GPUTexture[]
    uniform: Float32Array
    instances: number
    version: string
    textureVersions: string[]
}

/** Native 3D volume upload and ray-marching pipelines. */
export class WebGPUVolumeRenderer {
    readonly stats = { gridUploads: 0, transferUploads: 0, paletteUploads: 0 };
    private readonly items = new Map<number, VolumeItem>();
    private readonly layout: GPUBindGroupLayout;
    private readonly pipelines: GPURenderPipeline[];
    private readonly weightedPipelines: GPURenderPipeline[];
    private readonly sampler: GPUSampler;
    private readonly transparentPipeline: GPURenderPipeline;
    private readonly outlinePipeline: GPURenderPipeline;
    private readonly pickPipeline: GPURenderPipeline;
    private readonly colorDepthPipeline: GPURenderPipeline;
    private readonly peelingPipelines: Record<'near' | 'far' | 'front' | 'back', GPURenderPipeline>;

    constructor(private readonly context: WebGPUContext, shader: GPUShaderModule, cameraLayout: GPUBindGroupLayout, peeling: WebGPUDepthPeeling) {
        const { device } = context;
        this.layout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'read-only-storage' } },
            { binding: 2, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'unfilterable-float', viewDimension: '3d' } },
            { binding: 3, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 4, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 5, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
            { binding: 6, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'depth' } },
            ...[7, 8, 9, 10, 11, 12, 13].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'read-only-storage' as const } })),
        ] });
        this.pipelines = [true, false].map(emission => device.createRenderPipeline({
            label: 'molstar-volume-raymarch',
            layout: device.createPipelineLayout({ bindGroupLayouts: [cameraLayout, this.layout] }),
            vertex: { module: shader, entryPoint: 'vs' },
            fragment: { module: shader, entryPoint: 'fs', targets: [
                { format: context.format, blend: { color: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' }, alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' } } },
                { format: 'rgba32uint' },
                { format: 'rgba16float', writeMask: emission ? 15 : 8, blend: { color: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' }, alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' } } },
            ] },
            primitive: { topology: 'triangle-list' },
            depthStencil: { format: 'depth32float', depthWriteEnabled: false, depthCompare: 'less-equal' },
        }));
        this.weightedPipelines = [true, false].map(emission => device.createRenderPipeline({
            label: 'molstar-volume-weighted', layout: device.createPipelineLayout({ bindGroupLayouts: [cameraLayout, this.layout] }),
            vertex: { module: shader, entryPoint: 'vs' }, fragment: { module: shader, entryPoint: 'weighted', targets: WebGPUWeightedTransparency.targets(emission) },
            primitive: { topology: 'triangle-list' }, depthStencil: { format: 'depth32float', depthWriteEnabled: false, depthCompare: 'less-equal' },
        }));
        this.pickPipeline = device.createRenderPipeline({
            label: 'molstar-volume-picking', layout: device.createPipelineLayout({ bindGroupLayouts: [cameraLayout, this.layout] }),
            vertex: { module: shader, entryPoint: 'vs' }, fragment: { module: shader, entryPoint: 'pick', targets: [{ format: 'rgba32uint' }] },
            primitive: { topology: 'triangle-list' }, depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'less-equal' },
        });
        this.colorDepthPipeline = device.createRenderPipeline({
            label: 'molstar-volume-color-depth', layout: device.createPipelineLayout({ bindGroupLayouts: [cameraLayout, this.layout] }),
            vertex: { module: shader, entryPoint: 'vs' }, fragment: { module: shader, entryPoint: 'colorDepthId', targets: [{ format: 'rgba32uint' }] },
            primitive: { topology: 'triangle-list' }, depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'less-equal' },
        });
        const peelLayout = device.createPipelineLayout({ bindGroupLayouts: [cameraLayout, this.layout, peeling.emptyLayout, peeling.layout] });
        const peel = (phase: 'near' | 'far' | 'front' | 'back') => {
            const depth = phase === 'near' || phase === 'far';
            return device.createRenderPipeline({ label: `molstar-volume-peel-${phase}`, layout: peelLayout,
                vertex: { module: shader, entryPoint: 'vs' }, fragment: { module: shader, entryPoint: depth ? 'peelDepth' : phase === 'front' ? 'peelFront' : 'peelBack', targets: depth ? [] : WebGPUDepthPeeling.targets() },
                primitive: { topology: 'triangle-list' }, depthStencil: depth ? { format: 'depth32float', depthWriteEnabled: true, depthCompare: phase === 'near' ? 'less' : 'greater' } : undefined,
            });
        };
        this.peelingPipelines = { near: peel('near'), far: peel('far'), front: peel('front'), back: peel('back') };
        this.transparentPipeline = device.createRenderPipeline({ label: 'molstar-volume-transparent-color', layout: device.createPipelineLayout({ bindGroupLayouts: [cameraLayout, this.layout] }),
            vertex: { module: shader, entryPoint: 'vs' }, fragment: { module: shader, entryPoint: 'transparentColor', targets: [{ format: context.format, blend: { color: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' }, alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' } } }] },
            primitive: { topology: 'triangle-list' }, depthStencil: { format: 'depth32float', depthWriteEnabled: false, depthCompare: 'less-equal' },
        });
        this.outlinePipeline = device.createRenderPipeline({
            label: 'molstar-volume-outline-depth', layout: device.createPipelineLayout({ bindGroupLayouts: [cameraLayout, this.layout] }),
            vertex: { module: shader, entryPoint: 'vs' }, fragment: { module: shader, entryPoint: 'outlineDepth', targets: [{ format: 'rgba32float' }] },
            primitive: { topology: 'triangle-list' }, depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'less-equal' },
        });
        this.sampler = device.createSampler({ minFilter: 'linear', magFilter: 'linear' });
    }

    private version(object: GraphicsRenderObject) {
        const source = value<WebGPUTextureData | undefined>(object.values, 'tGridTex', undefined);
        return Object.values(object.values).map(cell => `${cell.ref.id}:${cell.ref.version}`).join(',') + `/${object.state.alphaFactor}/${source?.id}:${source?.version}`;
    }

    private createItem(object: GraphicsRenderObject, previous?: VolumeItem): VolumeItem {
        const { context } = this;
        const { device } = context;
        const v: RenderableValues = object.values;
        const source = value<WebGPUTextureData | undefined>(v, 'tGridTex', undefined);
        if (!(source instanceof WebGPUTextureData)) throw new Error('WebGPU volume rendering requires CPU-owned volume texture data.');
        const data = source.data;
        if (Math.max(data.width, data.height, data.depth) > device.limits.maxTextureDimension3D) throw new Error('Volume dimensions exceed the WebGPU device 3D texture limit.');
        const textureVersions = [`${v.tGridTex.ref.id}:${v.tGridTex.ref.version}/${source.id}/${source.version}`, `${v.tTransferTex.ref.id}:${v.tTransferTex.ref.version}`, `${v.tPalette.ref.id}:${v.tPalette.ref.version}`];
        const textures: GPUTexture[] = [], allocated: GPUTexture[] = [], buffers: GPUBuffer[] = [];
        const reuse = (index: number) => previous?.textureVersions[index] === textureVersions[index] ? previous.textures[index] : undefined;
        try {
            let grid = reuse(0);
            if (!grid) {
                grid = device.createTexture({ label: 'molstar-density-grid', size: { width: data.width, height: data.height, depthOrArrayLayers: data.depth }, dimension: '3d', format: 'rgba32float', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST });
                allocated.push(grid);
                const array = data.array instanceof Float32Array ? data.array : Float32Array.from(data.array, n => data.array instanceof Uint16Array ? fromHalfFloat(n) : n / 255);
                device.queue.writeTexture({ texture: grid }, array, { bytesPerRow: data.width * 16, rowsPerImage: data.height }, { width: data.width, height: data.height, depthOrArrayLayers: data.depth });
                this.stats.gridUploads++;
            }
            textures.push(grid);
            const transfer = value(v, 'tTransferTex', { array: new Uint8Array([0]), width: 1, height: 1 });
            let tf = reuse(1);
            if (!tf) {
                tf = device.createTexture({ size: { width: transfer.width, height: transfer.height }, format: 'r8unorm', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST });
                allocated.push(tf);
                device.queue.writeTexture({ texture: tf }, transfer.array, { bytesPerRow: transfer.width }, { width: transfer.width, height: transfer.height });
                this.stats.transferUploads++;
            }
            textures.push(tf);
            const palette = value(v, 'tPalette', { array: new Uint8Array([255, 255, 255]), width: 1, height: 1 });
            let pt = reuse(2);
            if (!pt) {
                const rgba = new Uint8Array(palette.width * palette.height * 4);
                for (let i = 0; i < rgba.length / 4; i++) { rgba.set(palette.array.subarray(i * 3, i * 3 + 3), i * 4); rgba[i * 4 + 3] = 255; }
                pt = device.createTexture({ size: { width: palette.width, height: palette.height }, format: 'rgba8unorm', usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST });
                allocated.push(pt);
                device.queue.writeTexture({ texture: pt }, rgba, { bytesPerRow: palette.width * 4 }, { width: palette.width, height: palette.height });
                this.stats.paletteUploads++;
            }
            textures.push(pt);

            const instances = value(v, 'instanceCount', 1);
            const transforms = value(v, 'aTransform', new Float32Array(16));
            const instanceIds = value(v, 'aInstance', new Float32Array(0));
            const unitToCartn = value(v, 'uUnitToCartn', Mat4.identity());
            const instanceData = new Float32Array(instances * 36);
            const instanceIntegers = new Uint32Array(instanceData.buffer);
            const world = Mat4(), inverse = Mat4();
            for (let i = 0; i < instances; i++) {
                Mat4.mul(world, Mat4.fromArray(Mat4(), transforms, i * 16), unitToCartn);
                if (!world.every(Number.isFinite) || !Mat4.tryInvert(inverse, world)) throw new Error('WebGPU volume instance transform is singular or non-finite.');
                instanceData.set(world, i * 36); instanceData.set(inverse, i * 36 + 16);
                instanceIntegers[i * 36 + 32] = instanceIds[i] ?? i;
            }
            const uniform = new Float32Array(60);
            const dimensions = value(v, 'uGridDim', [data.width, data.height, data.depth]);
            uniform.set(dimensions, 16);
            const order = value(v, 'dAxisOrder', '210').split('').map(Number);
            const strides = [0, 0, 0];
            strides[order[2]] = 1; strides[order[1]] = dimensions[order[2]]; strides[order[0]] = dimensions[order[2]] * dimensions[order[1]];
            uniform.set(strides, 20);
            uniform.set([value(v, 'uMaxSteps', 512), value(v, 'uTransferScale', 1 / 3), Math.max(0, Math.min(1, value(v, 'alpha', value(v, 'uAlpha', 1)) * object.state.alphaFactor)), value(v, 'dIgnoreLight', false) ? 1 : 0], 24);
            uniform.set(value(v, 'uColor', [1, 1, 1]), 28);
            uniform.set(value(v, 'uPaletteDomain', [0, 1]), 32);
            uniform[34] = value(v, 'tPalette', { filter: 'nearest' }).filter === 'linear' ? 1 : 0;
            const colorTypes: Record<string, number> = { uniform: 0, direct: 1, instance: 2, group: 3, groupInstance: 4, vertex: 5, vertexInstance: 6 };
            const colorType = colorTypes[value(v, 'dColorType', 'uniform')];
            if (colorType === undefined) throw new Error('WebGPU volume color granularity has not been implemented.');
            const ints = new Uint32Array(uniform.buffer);
            ints.set([object.id, value(v, 'uGroupCount', 1), colorType, value(v, 'uVertexCount', dimensions[0] * dimensions[1] * dimensions[2])], 36);
            uniform.set([value(v, 'uMarker', 0), value(v, 'dOverpaint', false) ? value(v, 'uOverpaintStrength', 1) : 0,
                value(v, 'dTransparency', false) ? value(v, 'uTransparencyStrength', 1) : 0,
                value(v, 'dEmissive', false) ? value(v, 'uEmissiveStrength', 1) : 0], 40);
            for (const [i, type] of ['dMarkerType', 'dOverpaintType', 'dTransparencyType', 'dEmissiveType'].entries()) {
                const granularity = value<string>(v, type, 'groupInstance');
                if (granularity !== 'instance' && granularity !== 'groupInstance' && granularity !== 'vertexInstance') throw new Error(`WebGPU volume ${type} '${granularity}' has not been implemented.`);
                ints[44 + i] = granularity === 'instance' ? 0 : granularity === 'vertexInstance' ? 2 : 1;
            }
            ints.set([value(v, 'dClipObjectCount', 0), value<string>(v, 'dClipVariant', 'pixel') === 'instance' ? 1 : 0, value<string>(v, 'dClippingType', 'groupInstance') === 'instance' ? 1 : 0, value(v, 'dClipping', false) ? 1 : 0], 48);
            uniform[52] = value(v, 'uEmissive', 0);
            uniform.set([value(v, 'uMetalness', 0), value(v, 'uRoughness', 1), value(v, 'uBumpiness', 0), value(v, 'dCelShaded', false) ? 1 : 0], 56);
            buffers.push(context.createBuffer(uniform, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST));
            buffers.push(context.createBuffer(instanceData, GPUBufferUsage.STORAGE));
            for (const name of ['tMarker', 'tTransparency', 'tOverpaint', 'tEmissive', 'tColor']) {
                const bytes = value(v, name, { array: new Uint8Array(0) }).array;
                if (bytes.byteLength > device.limits.maxStorageBufferBindingSize) throw new Error('WebGPU volume theme exceeds the device storage buffer binding limit.');
                buffers.push(context.createBuffer(bytes, GPUBufferUsage.STORAGE));
            }
            buffers.push(context.createBuffer(createClipData(v), GPUBufferUsage.STORAGE));
            buffers.push(context.createBuffer(value(v, 'tClipping', { array: new Uint8Array(1) }).array, GPUBufferUsage.STORAGE));
            return { object, buffers, textures, uniform, instances, version: this.version(object), textureVersions };
        } catch (error) {
            for (const buffer of buffers) buffer.destroy();
            for (const texture of allocated) texture.destroy();
            throw error;
        }
    }

    sync(objects: readonly GraphicsRenderObject[]) {
        const ids = new Set(objects.map(object => object.id));
        for (const [id, item] of this.items) if (!ids.has(id)) { this.destroyItem(item); this.items.delete(id); }
    }

    render(pass: GPURenderPassEncoder, objects: readonly GraphicsRenderObject[], camera: Camera, cameraBindGroup: GPUBindGroup, depth: GPUTexture, emission = true, picking = false, outline = false, transparent = false, weighted = false, colorDepth = false, peeling?: 'near' | 'far' | 'front' | 'back'): number {
        const inverseCamera = Mat4.invert(Mat4(), Mat4.mul(Mat4(), camera.projection, getWebGPUModelView(camera)));
        pass.setPipeline(peeling ? this.peelingPipelines[peeling] : colorDepth ? this.colorDepthPipeline : weighted ? this.weightedPipelines[emission ? 0 : 1] : outline ? this.outlinePipeline : transparent ? this.transparentPipeline : picking ? this.pickPipeline : this.pipelines[emission ? 0 : 1]); pass.setBindGroup(0, cameraBindGroup);
        let count = 0;
        for (const object of objects) {
            if (picking && (!object.state.pickable || object.state.colorOnly)) continue;
            if (!object.state.visible || !value(object.values, 'instanceCount', 0)) continue;
            let item = this.items.get(object.id);
            if (!item || item.version !== this.version(object)) {
                const next = this.createItem(object, item);
                if (item) this.destroyItem(item, next.textures);
                this.items.set(object.id, next); item = next;
            }
            item.uniform.set(inverseCamera, 0);
            this.context.device.queue.writeBuffer(item.buffers[0], 0, item.uniform);
            const bindGroup = this.context.device.createBindGroup({ layout: this.layout, entries: [
                { binding: 0, resource: { buffer: item.buffers[0] } }, { binding: 1, resource: { buffer: item.buffers[1] } },
                ...[2, 3, 4].map((binding, i) => ({ binding, resource: item!.textures[i].createView() })),
                { binding: 5, resource: this.sampler }, { binding: 6, resource: depth.createView() },
                ...[7, 8, 9, 10, 11, 12, 13].map((binding, i) => ({ binding, resource: { buffer: item!.buffers[i + 2] } })),
            ] });
            pass.setBindGroup(1, bindGroup); pass.draw(3, item.instances);
            count++;
        }
        return count;
    }

    getObject(id: number) { return this.items.get(id)?.object; }
    private destroyItem(item: VolumeItem, retained: readonly GPUTexture[] = []) {
        for (const buffer of item.buffers) buffer.destroy();
        for (const texture of item.textures) if (!retained.includes(texture)) texture.destroy();
    }
    dispose() { for (const item of this.items.values()) this.destroyItem(item); this.items.clear(); }
}
