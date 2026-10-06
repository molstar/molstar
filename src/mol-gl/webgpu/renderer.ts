/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { IlluminationProps } from '../../mol-canvas3d/passes/illumination';
import { createWebGPUCameraHelper } from '../../mol-canvas3d/helper/camera-helper';
import { createWebGPUHandleHelper } from '../../mol-canvas3d/helper/handle-helper';
import { createWebGPUPointerHelper } from '../../mol-canvas3d/helper/pointer-helper';
import { StereoCameraProps } from '../../mol-canvas3d/camera/stereo';
import { WebGPUStereoCamera } from './stereo';
import { WebGPUTracingQuality } from './tracing-quality';
import { WebGPUTracing } from './tracing';
import { WebGPUIlluminationCompose } from './illumination-compose';
import { WebGPUTracingInput } from './tracing-input';
import { createLightingData, WebGPULightingByteSize } from './lighting';
import { smaaShader, WebGPUSmaa } from './smaa';
import { bloomShader } from './bloom';
import { WebGPUPostprocessing } from './postprocessing';
import { WebGPUMarking } from './marking';
import { WebGPUBackground } from './background';
import { WebGPUMultiSample } from './multi-sample';
import { MultiSampleProps } from '../../mol-canvas3d/passes/multi-sample';
import { AssetManager } from '../../mol-util/assets';
import { MarkingProps } from '../../mol-canvas3d/passes/marking';
import { postprocessingShader } from './postprocessing-shader';
import { PostprocessingProps } from '../../mol-canvas3d/passes/postprocessing';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { getWebGPUModelView } from './camera';
import { Camera } from '../../mol-canvas3d/camera';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Color } from '../../mol-util/color';
import { GraphicsRenderObject, isMergedRenderObject } from '../render-object';
import { WebGPUDebugRegistry } from '../../mol-canvas3d/helper/debug-registry';
import { RendererProps } from '../renderer';
import { RenderableValues } from '../renderable/schema';
import { WebGPUContext } from './context';
import { createWebGPUGeometry, createWebGPUThemes, value, WebGPUGeometry, WebGPUThemeStride, createWebGPUOverpaintThemes } from './geometry';
import { geometryShader, rasterGeometryShader } from './shader';
import { GPUBufferUsage, GPUShaderStage, GPUTextureUsage } from './compat';
import { WebGPUVolumeRenderer } from './volume';
import { volumeShader } from './volume-shader';
import { createImageStyles } from './image';
import { createClipData } from './clip';
import { WebGPUTextureData } from './texture-data';
import { WebGPUWeightedTransparency } from './transparency';
import { WebGPUDepthPeeling } from './depth-peeling';
import { getSphereLodInstanceRanges, SphereInstanceRange } from './sphere-lod';

interface RenderItem {
    object: GraphicsRenderObject
    buffers: GPUBuffer[]
    borrowedGeometry?: GPUBuffer
    texture: GPUTexture
    imageTextures: GPUTexture[]
    ownedImageTextures: GPUTexture[]
    bindGroup: GPUBindGroup
    count: number
    instances: number
    version: string
    transparent: boolean
    surface: boolean
    hasOpaque: boolean
    geometry: WebGPUGeometry
    geometryVersion: string
    transformVersion: string
    textureVersions: string[]
    themeVersion: string
    themeHasTransparency: boolean
    themeHasOpaque: boolean
}

export interface WebGPUPick {
    id: PickingId
    /** WebGPU normalized depth in 0..1. */
    depth: number
}

/** Native WGSL rendering of Mol* render objects. Coordinates for readback are top-left. */
export class WebGPURenderer {
    private readonly items = new Map<number, RenderItem>();
    private readonly pipelines = new Map<string, GPURenderPipeline>();
    private readonly cameraBuffer: GPUBuffer;
    private readonly lightingBuffer: GPUBuffer;
    private readonly cameraBindGroup: GPUBindGroup;
    private readonly helperCameraBuffer: GPUBuffer;
    private readonly helperCameraBindGroup: GPUBindGroup;
    private readonly pointerCameraBuffer: GPUBuffer;
    private readonly pointerCameraBindGroup: GPUBindGroup;
    private helperDepth?: GPUTexture;
    private cameraHelper?: ReturnType<typeof createWebGPUCameraHelper>;
    private readonly helperIds = new Set<number>();
    setCameraHelper(helper: ReturnType<typeof createWebGPUCameraHelper> | undefined) { this.cameraHelper = helper; }
    private handleHelper?: ReturnType<typeof createWebGPUHandleHelper>;
    private pointerHelper?: ReturnType<typeof createWebGPUPointerHelper>;
    setWorldHelpers(handle: ReturnType<typeof createWebGPUHandleHelper> | undefined, pointer: ReturnType<typeof createWebGPUPointerHelper> | undefined) { this.handleHelper = handle; this.pointerHelper = pointer; }
    private debugRegistry?: WebGPUDebugRegistry;
    setDebugRegistry(registry: WebGPUDebugRegistry | undefined) { this.debugRegistry = registry; }
    private readonly objectLayout: GPUBindGroupLayout;
    private readonly pipelineLayout: GPUPipelineLayout;
    private readonly shader: GPUShaderModule;
    private readonly rasterShader: GPUShaderModule;
    private readonly volumes: WebGPUVolumeRenderer;
    private readonly postprocessing: WebGPUPostprocessing;
    private readonly marking: WebGPUMarking;
    private readonly background: WebGPUBackground;
    private readonly weightedTransparency: WebGPUWeightedTransparency;
    private readonly depthPeeling: WebGPUDepthPeeling;
    private readonly peelingPipelineLayout: GPUPipelineLayout;
    private transparency: 'blended' | 'wboit' | 'dpoit' = 'blended';
    private peelingIterations = 2;
    setDpoitIterations(iterations: number) {
        const value = Math.max(1, Math.min(10, Math.round(iterations)));
        if (!Number.isFinite(value)) throw new Error('Depth-peeling iterations must be finite.');
        if (value === this.peelingIterations) return;
        this.peelingIterations = value; this.multiSample.reset(); this.illuminationKey = '';
        for (const sample of this.stereo?.samples ?? []) sample.reset();
    }
    get transparencyMode() { return this.transparency; }
    setTransparency(mode: 'blended' | 'wboit' | 'dpoit') {
        if (mode === this.transparency) return;
        this.transparency = mode; this.multiSample.reset(); this.illuminationKey = '';
        for (const sample of this.stereo?.samples ?? []) sample.reset();
    }
    private multiSample: WebGPUMultiSample;
    private readonly monoMultiSample: WebGPUMultiSample;
    private stereo?: { camera: WebGPUStereoCamera, samples: [WebGPUMultiSample, WebGPUMultiSample], color?: GPUTexture, picking?: GPUTexture, key?: string };
    private stereoActive = false;
    private outputPicking?: GPUTexture;
    private readonly tracingInput: WebGPUTracingInput;
    private illuminationIteration = 0;
    private illuminationMax = 0;
    private illuminationKey = '';
    private readonly tracingQuality = new WebGPUTracingQuality();
    private quality = { rendersPerFrame: 1, steps: 1, refineSteps: 0 };
    get illuminationQuality() { return { ...this.quality }; }
    get illuminationNeedsFrame() { return this.illuminationIteration < this.illuminationMax; }
    get illuminationProgress() { return this.illuminationIteration; }
    get multiSampleNeedsFrame() { return this.stereoActive ? this.stereo!.samples.some(s => s.needsFrame) : this.multiSample.needsFrame; }
    getPickingCamera(x: number, y: number, fallback: Camera): Camera {
        if (!this.stereoActive) return fallback;
        const { left, right } = this.stereo!.camera;
        for (const eye of [left, right]) {
            const v = eye.viewport;
            if (x >= v.x && x < v.x + v.width && y >= this.height - v.y - v.height && y < this.height - v.y) return eye;
        }
        return fallback;
    }
    get backgroundChanged() { return this.background.changed; }
    async updateBackground(props?: PostprocessingProps) { await this.background.ready(props?.enabled ? props.background : { variant: { name: 'off', params: {} } }); }
    private readonly markingPipelineLayout: GPUPipelineLayout;
    private outputColor?: GPUTexture;
    private color?: GPUTexture;
    private emissive?: GPUTexture;
    private transparentColor?: GPUTexture;
    private outlineDepthAlpha?: GPUTexture;
    private outlineDepth?: GPUTexture;
    private depth?: GPUTexture;
    private colorPickDepth?: GPUTexture;
    private opaqueDepth?: GPUTexture;
    private pickTexture?: GPUTexture;
    private selectionTexture?: GPUTexture;
    private selectionDepth?: GPUTexture;
    private selectionOpaqueDepth?: GPUTexture;
    private width = 0;
    private height = 0;
    private disposed = false;
    private time = 0;
    get animationTime() { return this.time; }
    get volumeResourceStats() { return { ...this.volumes.stats }; }
    readonly geometryResourceStats = { builds: 0, vertexUploads: 0, indexUploads: 0, transformUploads: 0 };
    readonly textureResourceStats = { uploads: 0 };
    readonly bufferResourceStats = { allocations: 0, updates: 0, uploadedBytes: 0 };
    readonly themeResourceStats = { builds: 0, uploads: 0 };
    readonly sphereLodStats = { draws: 0, indices: 0 };
    private readonly sphereLodRanges = new Map<RenderItem, SphereInstanceRange[][]>();
    private sphereAnimation = false;
    setTime(time: number) { this.time = time; }
    readonly stats = { drawCount: 0, instanceCount: 0, triangleCount: 0 };

    private constructor(readonly context: WebGPUContext, shader: GPUShaderModule, volumeModule: GPUShaderModule, postModule: GPUShaderModule, bloomModule: GPUShaderModule, smaa: WebGPUSmaa, private readonly tracing: WebGPUTracing, private readonly illuminationCompose: WebGPUIlluminationCompose, assets?: AssetManager) {
        const { device } = context;
        this.shader = shader;
        this.rasterShader = device.createShaderModule({ label: 'molstar-native-raster-geometry', code: rasterGeometryShader });
        const cameraLayout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
        ] });
        this.objectLayout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            ...[1, 2, 3].map(binding => ({ binding, visibility: GPUShaderStage.VERTEX | (binding === 2 ? GPUShaderStage.FRAGMENT : 0), buffer: { type: 'read-only-storage' as const } })),
            { binding: 4, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 5, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
            { binding: 6, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 7, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'unfilterable-float' } },
            { binding: 8, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 9, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'read-only-storage' } },
            ...[10, 11, 12].map(binding => ({ binding, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'read-only-storage' as const } })),
            ...[13, 14, 15, 16, 17].map(binding => ({ binding, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' as const } })),
            { binding: 18, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'read-only-storage' } },
        ] });
        this.pipelineLayout = device.createPipelineLayout({ bindGroupLayouts: [cameraLayout, this.objectLayout] });
        this.marking = new WebGPUMarking(context);
        this.background = new WebGPUBackground(context, assets);
        this.weightedTransparency = new WebGPUWeightedTransparency(context);
        this.depthPeeling = new WebGPUDepthPeeling(context);
        this.peelingPipelineLayout = device.createPipelineLayout({ bindGroupLayouts: [cameraLayout, this.objectLayout, this.depthPeeling.emptyLayout, this.depthPeeling.layout] });
        this.multiSample = new WebGPUMultiSample(context);
        this.monoMultiSample = this.multiSample;
        this.tracingInput = new WebGPUTracingInput(context, this.rasterShader, this.pipelineLayout, shader);
        this.markingPipelineLayout = device.createPipelineLayout({ bindGroupLayouts: [cameraLayout, this.objectLayout, this.marking.geometryLayout] });
        this.cameraBuffer = device.createBuffer({ size: 256, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.lightingBuffer = device.createBuffer({ size: WebGPULightingByteSize, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.cameraBindGroup = device.createBindGroup({ layout: cameraLayout, entries: [{ binding: 0, resource: { buffer: this.cameraBuffer } }, { binding: 1, resource: { buffer: this.lightingBuffer } }] });
        this.helperCameraBuffer = device.createBuffer({ size: 256, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.helperCameraBindGroup = device.createBindGroup({ layout: cameraLayout, entries: [{ binding: 0, resource: { buffer: this.helperCameraBuffer } }, { binding: 1, resource: { buffer: this.lightingBuffer } }] });
        this.pointerCameraBuffer = device.createBuffer({ size: 256, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
        this.pointerCameraBindGroup = device.createBindGroup({ layout: cameraLayout, entries: [{ binding: 0, resource: { buffer: this.pointerCameraBuffer } }, { binding: 1, resource: { buffer: this.lightingBuffer } }] });
        this.volumes = new WebGPUVolumeRenderer(context, volumeModule, cameraLayout, this.depthPeeling);
        this.postprocessing = new WebGPUPostprocessing(context, postModule, bloomModule, smaa);
    }

    static async create(context: WebGPUContext, assets?: AssetManager) {
        const shader = context.device.createShaderModule({ label: 'molstar-native-geometry', code: geometryShader });
        const volumeModule = context.device.createShaderModule({ label: 'molstar-native-volume', code: volumeShader });
        const postModule = context.device.createShaderModule({ label: 'molstar-postprocessing', code: postprocessingShader });
        const bloomModule = context.device.createShaderModule({ label: 'molstar-bloom', code: bloomShader });
        const smaaModule = context.device.createShaderModule({ label: 'molstar-smaa', code: smaaShader });
        const infos = await Promise.all([shader.getCompilationInfo(), volumeModule.getCompilationInfo(), postModule.getCompilationInfo(), bloomModule.getCompilationInfo(), smaaModule.getCompilationInfo()]);
        const errors = infos.flatMap(info => info.messages).filter(m => m.type === 'error');
        if (errors.length) throw new Error(`WebGPU shader compilation failed:\n${errors.map(e => `${e.lineNum}:${e.linePos}: ${e.message}`).join('\n')}`);
        const [smaaPass, tracing, compose] = await Promise.all([WebGPUSmaa.create(context, smaaModule), WebGPUTracing.create(context), WebGPUIlluminationCompose.create(context)]);
        return new WebGPURenderer(context, shader, volumeModule, postModule, bloomModule, smaaPass, tracing, compose, assets);
    }

    private pipeline(transparent: boolean, cull: boolean, emission = true, flipSided = false, phase: 'all' | 'opaque' | 'transparent' = 'all', sphere = false): GPURenderPipeline {
        const key = `${transparent}/${cull}/${emission}/${flipSided}/${phase}/${sphere}`;
        const module = sphere ? this.shader : this.rasterShader;
        let pipeline = this.pipelines.get(key);
        if (!pipeline) {
            pipeline = this.context.device.createRenderPipeline({
                label: `molstar-geometry-${key}`,
                layout: this.pipelineLayout,
                vertex: { module, entryPoint: 'vs' },
                fragment: { module, entryPoint: phase === 'opaque' ? 'fsOpaque' : phase === 'transparent' ? 'fsTransparent' : 'fs', targets: [
                    { format: this.context.format, blend: transparent ? {
                        color: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha', operation: 'add' },
                        alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha', operation: 'add' },
                    } : undefined },
                    { format: 'rgba32uint' },
                    { format: 'rgba16float', writeMask: emission ? 15 : 8, blend: transparent ? { color: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' }, alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' } } : undefined },
                ] },
                primitive: { topology: 'triangle-list', cullMode: cull ? (flipSided ? 'front' : 'back') : 'none' },
                depthStencil: { format: 'depth32float', depthWriteEnabled: !transparent, depthCompare: 'less-equal' },
            });
            this.pipelines.set(key, pipeline);
        }
        return pipeline;
    }

    private pickingPipeline(cull: boolean, flipSided: boolean, sphere = false) {
        const key = `pick/${cull}/${flipSided}/${sphere}`;
        const module = sphere ? this.shader : this.rasterShader;
        let pipeline = this.pipelines.get(key);
        if (!pipeline) {
            pipeline = this.context.device.createRenderPipeline({ label: `molstar-${key}`, layout: this.pipelineLayout,
                vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'pick', targets: [{ format: 'rgba32uint' }] },
                primitive: { topology: 'triangle-list', cullMode: cull ? (flipSided ? 'front' : 'back') : 'none' },
                depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'less-equal' },
            });
            this.pipelines.set(key, pipeline);
        }
        return pipeline;
    }
    private weightedPipeline(cull: boolean, flipSided: boolean, emission: boolean, sphere = false) {
        const key = `weighted/${cull}/${flipSided}/${emission}/${sphere}`;
        const module = sphere ? this.shader : this.rasterShader;
        let pipeline = this.pipelines.get(key);
        if (!pipeline) {
            pipeline = this.context.device.createRenderPipeline({ label: `molstar-${key}`, layout: this.pipelineLayout,
                vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'weighted', targets: WebGPUWeightedTransparency.targets(emission) },
                primitive: { topology: 'triangle-list', cullMode: cull ? (flipSided ? 'front' : 'back') : 'none' },
                depthStencil: { format: 'depth32float', depthWriteEnabled: false, depthCompare: 'less-equal' },
            });
            this.pipelines.set(key, pipeline);
        }
        return pipeline;
    }
    private colorDepthPipeline(cull: boolean, flipSided: boolean, sphere = false) {
        const key = `color-depth/${cull}/${flipSided}/${sphere}`;
        const module = sphere ? this.shader : this.rasterShader;
        let pipeline = this.pipelines.get(key);
        if (!pipeline) {
            pipeline = this.context.device.createRenderPipeline({ label: `molstar-${key}`, layout: this.pipelineLayout,
                vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: 'colorDepthId', targets: [{ format: 'rgba32uint' }] },
                primitive: { topology: 'triangle-list', cullMode: cull ? (flipSided ? 'front' : 'back') : 'none' },
                depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'less-equal' },
            }); this.pipelines.set(key, pipeline);
        } return pipeline;
    }
    private peelingPipeline(phase: 'near' | 'far' | 'front' | 'back', cull: boolean, flipSided: boolean, sphere = false) {
        const key = `peel/${phase}/${cull}/${flipSided}/${sphere}`;
        const module = sphere ? this.shader : this.rasterShader;
        let pipeline = this.pipelines.get(key);
        if (!pipeline) {
            const depth = phase === 'near' || phase === 'far';
            pipeline = this.context.device.createRenderPipeline({ label: `molstar-${key}`, layout: this.peelingPipelineLayout,
                vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: depth ? 'peelDepth' : phase === 'front' ? 'peelFront' : 'peelBack', targets: depth ? [] : WebGPUDepthPeeling.targets() },
                primitive: { topology: 'triangle-list', cullMode: cull ? (flipSided ? 'front' : 'back') : 'none' },
                depthStencil: depth ? { format: 'depth32float', depthWriteEnabled: true, depthCompare: phase === 'near' ? 'less' : 'greater' } : undefined,
            }); this.pipelines.set(key, pipeline);
        } return pipeline;
    }

    private markingPipeline(mask: boolean, cull: boolean, flipSided: boolean, sphere = false) {
        const key = `marking/${mask}/${cull}/${flipSided}/${sphere}`;
        const module = sphere ? this.shader : this.rasterShader;
        let pipeline = this.pipelines.get(key);
        if (!pipeline) {
            pipeline = this.context.device.createRenderPipeline({ label: `molstar-${key}`, layout: mask ? this.markingPipelineLayout : this.pipelineLayout,
                vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: mask ? 'markingMask' : 'markingDepth', targets: mask ? [{ format: 'rgba16float' }] : [] },
                primitive: { topology: 'triangle-list', cullMode: cull ? (flipSided ? 'front' : 'back') : 'none' },
                depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'less-equal' },
            }); this.pipelines.set(key, pipeline);
        }
        return pipeline;
    }

    private transparentColorPipeline(cull: boolean, flipSided: boolean, surface: boolean, sphere = false) {
        const key = `transparent-color/${cull}/${flipSided}/${surface}/${sphere}`;
        const module = sphere ? this.shader : this.rasterShader;
        let pipeline = this.pipelines.get(key);
        if (!pipeline) {
            pipeline = this.context.device.createRenderPipeline({ label: `molstar-${key}`, layout: this.pipelineLayout,
                vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: surface ? 'transparentSurfaceColor' : 'transparentColor', targets: [{ format: this.context.format, blend: { color: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' }, alpha: { srcFactor: 'one', dstFactor: 'one-minus-src-alpha' } } }] },
                primitive: { topology: 'triangle-list', cullMode: cull ? (flipSided ? 'front' : 'back') : 'none' },
                depthStencil: { format: 'depth32float', depthWriteEnabled: false, depthCompare: 'less-equal' },
            }); this.pipelines.set(key, pipeline);
        }
        return pipeline;
    }

    private outlineDepthPipeline(cull: boolean, flipSided: boolean, surface: boolean, sphere = false) {
        const key = `outline-depth/${cull}/${flipSided}/${surface}/${sphere}`;
        const module = sphere ? this.shader : this.rasterShader;
        let pipeline = this.pipelines.get(key);
        if (!pipeline) {
            pipeline = this.context.device.createRenderPipeline({ label: `molstar-${key}`, layout: this.pipelineLayout,
                vertex: { module, entryPoint: 'vs' }, fragment: { module, entryPoint: surface ? 'outlineSurfaceDepth' : 'outlineDepth', targets: [{ format: 'rgba32float' }] },
                primitive: { topology: 'triangle-list', cullMode: cull ? (flipSided ? 'front' : 'back') : 'none' },
                depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'less-equal' },
            }); this.pipelines.set(key, pipeline);
        }
        return pipeline;
    }

    private version(object: GraphicsRenderObject) {
        return this.valueVersions(object.values, Object.keys(object.values)) + `/${object.state.alphaFactor}/${object.state.opaque}`;
    }

    private countDraw(indices: number, instances: number) {
        this.stats.drawCount++; this.stats.instanceCount += instances; this.stats.triangleCount += indices / 3 * instances;
    }

    private drawItem(pass: GPURenderPassEncoder, item: RenderItem, camera: Camera, countStats = false) {
        if (!item.geometry.sphereLods) { pass.drawIndexed(item.count, item.instances); if (countStats) this.countDraw(item.count, item.instances); return; }
        let ranges = this.sphereLodRanges.get(item);
        if (!ranges) {
            ranges = item.geometry.sphereLods.map(l => getSphereLodInstanceRanges(item.object.values, item.instances, camera, l.min, l.max, this.sphereAnimation));
            this.sphereLodRanges.set(item, ranges);
        }
        for (let i = 0; i < item.geometry.sphereLods.length; i++) {
            const level = item.geometry.sphereLods[i];
            if (!level.count) continue;
            for (const range of ranges[i]) {
                pass.drawIndexed(level.count, range.count, level.first, 0, range.first);
                if (countStats) this.countDraw(level.count, range.count);
                this.sphereLodStats.draws++; this.sphereLodStats.indices += level.count * range.count;
            }
        }
    }

    private valueVersions(values: RenderableValues, keys: readonly string[]) {
        return keys.map(key => {
            const cell = values[key];
            return cell ? `${key}:${cell.ref.id}:${cell.ref.version}:${cell.ref.value instanceof WebGPUTextureData ? `${cell.ref.value.id}:${cell.ref.value.version}` : ''}` : `${key}:-`;
        }).join(',');
    }

    private createItem(object: GraphicsRenderObject, previous?: RenderItem): RenderItem {
        const { context } = this;
        const { device } = context;
        const v: RenderableValues = object.values;
        // Geometry construction reads attributes and source textures, independently
        // of the material/marker/theme values expanded into the theme buffer.
        const geometryVersion = object.type + '/' + this.valueVersions(v, [
            'aPosition', 'aStart', 'aEnd', 'aNormal', 'aGroup', 'aMapping', 'aUv', 'aTexCoord', 'aDepth', 'elements',
            'centerBuffer', 'groupBuffer', 'aScale', 'aCap', 'aColorMode', 'drawCount', 'dSolidInterior',
            'tPosition', 'tNormal', 'tGroup', 'uGeoTexDim', 'meta',
        ]) + (object.type === 'spheres' && value<unknown[]>(v, 'lodLevels', []).length ? '/' + this.valueVersions(v, ['tPositionGroup', 'lodLevels']) : '');
        const reuseGeometry = previous?.geometryVersion === geometryVersion;
        const geometry = reuseGeometry ? previous.geometry : createWebGPUGeometry(object);
        if (!reuseGeometry) this.geometryResourceStats.builds++;
        const themeVersion = geometryVersion + '/' + this.valueVersions(v, object.type === 'image' ? ['instanceCount', 'aInstance', 'uAlpha'] : [
            'instanceCount', 'aInstance', 'uGroupCount', 'uVertexCount', 'tColor', 'dColorType', 'uColor', 'dDualColor',
            'tSize', 'dSizeType', 'uSize', 'uSizeFactor', 'tMarker', 'dMarkerType', 'uMarker', 'alpha', 'uAlpha',
            'tTransparency', 'dTransparency', 'dTransparencyType', 'uTransparencyStrength',
            'tEmissive', 'dEmissive', 'dEmissiveType', 'uEmissive', 'uEmissiveStrength',
            'tSubstance', 'dSubstance', 'dSubstanceType', 'uSubstanceStrength', 'uMetalness', 'uRoughness', 'uBumpiness',
        ]) + `/${object.state.alphaFactor}`;
        const reuseThemes = previous?.themeVersion === themeVersion;
        const themes = reuseThemes ? undefined : object.type === 'image' ? new Float32Array(value(v, 'instanceCount', 1) * geometry.vertexCount * WebGPUThemeStride) : createWebGPUThemes(object, geometry);
        if (!reuseThemes) this.themeResourceStats.builds++;
        if (object.type === 'image' && !reuseThemes) {
            const ids = value(v, 'aInstance', new Float32Array(0));
            for (let i = 0; i < value(v, 'instanceCount', 1); i++) for (let j = 0; j < geometry.vertexCount; j++) {
                const o = (i * geometry.vertexCount + j) * WebGPUThemeStride;
                themes![o + 3] = value(v, 'uAlpha', 1) * object.state.alphaFactor;
                themes![o + 7] = ids[i] ?? i;
            }
        }
        const instances = value(v, 'instanceCount', 1);
        const transforms = value(v, 'aTransform', new Float32Array(instances * 16));
        const transformVersion = this.valueVersions(v, ['aTransform', 'instanceCount']);
        const buffers: GPUBuffer[] = [], allocatedBuffers: GPUBuffer[] = [];
        const pendingWrites: { buffer: GPUBuffer, data: ArrayBufferView }[] = [];
        const allocateBuffer = (data: ArrayBufferView, usage: number) => {
            const buffer = context.createBuffer(data, usage);
            this.bufferResourceStats.allocations++;
            this.bufferResourceStats.uploadedBytes += buffer.size;
            allocatedBuffers.push(buffer); return buffer;
        };
        const updateBuffer = (index: number, data: ArrayBufferView, usage: number) => {
            const prior = previous?.buffers[index];
            if (prior && prior !== previous?.borrowedGeometry && (prior.usage & GPUBufferUsage.COPY_DST) && prior.size >= Math.max(4, data.byteLength)) {
                // Clear unused capacity when shrinking, including packed byte
                // channels whose shader bounds use the storage array length.
                let upload = data;
                if (data.byteLength !== prior.size) {
                    const padded = new Uint8Array(prior.size);
                    padded.set(new Uint8Array(data.buffer, data.byteOffset, data.byteLength)); upload = padded;
                }
                pendingWrites.push({ buffer: prior, data: upload }); return prior;
            }
            return allocateBuffer(data, usage | GPUBufferUsage.COPY_DST);
        };
        let texture: GPUTexture | undefined;
        const imageTextures: GPUTexture[] = [], ownedImageTextures: GPUTexture[] = [];
        const textureVersions: string[] = [], allocatedTextures: GPUTexture[] = [];
        const uploadTexture = (index: number, version: string, data: () => { width: number, height: number, array: ArrayBufferView }, format: GPUTextureFormat) => {
            textureVersions[index] = version;
            if (previous?.textureVersions[index] === version) return index === 0 ? previous.texture : previous.imageTextures[index - 1];
            const image = data();
            const result = device.createTexture({ size: { width: image.width, height: image.height }, format, usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST });
            allocatedTextures.push(result);
            device.queue.writeTexture({ texture: result }, image.array, { bytesPerRow: image.width * (format === 'r8unorm' ? 1 : 4) }, { width: image.width, height: image.height });
            this.textureResourceStats.uploads++;
            return result;
        };
        try {
            if (geometry.native && geometry.native.device !== device) throw new Error('Native geometry belongs to another GPU device.');
            for (const [index, data] of [geometry.vertices, transforms.subarray(0, instances * 16), themes ?? new Float32Array(0)].entries()) {
                if (data.byteLength > device.limits.maxStorageBufferBindingSize) throw new Error(`WebGPU render object ${object.id} exceeds the storage buffer binding limit.`);
                const reused = index === 0 && reuseGeometry ? previous.buffers[0]
                    : index === 1 && previous?.transformVersion === transformVersion ? previous.buffers[1]
                        : index === 2 && reuseThemes ? previous.buffers[2] : undefined;
                buffers.push(reused ?? (index === 0 && geometry.native ? geometry.native.buffer : updateBuffer(index, data, GPUBufferUsage.STORAGE)));
                if (!reused && index === 0 && !geometry.native) this.geometryResourceStats.vertexUploads++;
                if (!reused && index === 1) this.geometryResourceStats.transformUploads++;
            }
            buffers.push(reuseGeometry ? previous.buffers[3] : allocateBuffer(geometry.indices, GPUBufferUsage.INDEX));
            if (!reuseGeometry) this.geometryResourceStats.indexUploads++;
            const uniform = new Float32Array(144);
            new Uint32Array(uniform.buffer).set([object.id, geometry.vertexCount, geometry.kind, value(v, 'dClipObjectCount', 0)]);
            uniform.set([value(v, 'dIgnoreLight', false) ? 1 : 0, value(v, geometry.kind === 1 ? 'dLineSizeAttenuation' : 'dPointSizeAttenuation', false) ? 1 : 0,
                { square: 0, circle: 1, fuzzy: 2 }[value(v, 'dPointStyle', 'square') as 'square' | 'circle' | 'fuzzy'], value(v, 'dFlipSided', false) ? 1 : 0], 4);
            uniform.set([value(v, 'uOffsetX', 0), value(v, 'uOffsetY', 0), value(v, 'uOffsetZ', 0), 0], 8);
            uniform.set(value(v, 'uBorderColor', [0, 0, 0]), 12);
            uniform[15] = value(v, 'uBorderWidth', 0);
            uniform.set(value(v, 'uBackgroundColor', [0, 0, 0]), 16);
            uniform[19] = value(v, 'uBackgroundOpacity', 0);
            uniform.set([value(v, 'uIsoLevel', -1), { nearest: 0, catmulrom: 1, mitchell: 2, bspline: 3 }[value(v, 'dInterpolation', 'nearest') as 'nearest'], value(v, 'dUsePalette', false) ? 1 : 0, value(v, 'uGroupCount', 1)], 20);
            uniform.set(value(v, 'uPaletteDefault', [1, 1, 1]), 24);
            uniform.set(value(v, 'uTrimCenter', [0, 0, 0]), 28); uniform[31] = value(v, 'uTrimType', 0);
            uniform.set(value(v, 'uTrimRotation', [0, 0, 0, 1]), 32);
            uniform.set(value(v, 'uTrimScale', [1, 1, 1]), 36);
            uniform.set(value(v, 'uTrimTransform', Mat4.identity()), 40);
            uniform.set([value<string>(v, 'dClipVariant', 'pixel') === 'instance' ? 1 : 0, value<string>(v, 'dClippingType', 'groupInstance') === 'instance' ? 1 : 0, value(v, 'dClipping', false) ? 1 : 0, value(v, 'uGroupCount', 1)], 56);
            uniform.set(value(v, 'uInvariantBoundingSphere', [0, 0, 0, 0]), 60);
            uniform.set([value(v, 'uWiggleSpeed', 0), value(v, 'uWiggleAmplitude', 0), value(v, 'uWiggleFrequency', 0), value(v, 'uWiggleMode', 0)], 64);
            uniform.set([value(v, 'uTumbleSpeed', 0), value(v, 'uTumbleAmplitude', 0), value(v, 'uTumbleFrequency', 0), 0], 68);
            uniform.set([value(v, 'dWiggle', false) ? 1 : 0, value<string>(v, 'dWiggleType', 'groupInstance') === 'instance' ? 1 : 0, value(v, 'uWiggleStrength', 1), 0], 72);
            uniform.set([value(v, 'uMetalness', 0), value(v, 'uRoughness', 1), value(v, 'uBumpiness', 0), value(v, 'dCelShaded', false) ? 1 : 0], 76);
            uniform.set([value(v, 'uBumpFrequency', 0), value(v, 'uBumpAmplitude', 0), value(v, 'dFlatShaded', false) ? 1 : 0, { off: 0, on: 1, inverted: 2 }[value<string>(v, 'dXrayShaded', 'off') as 'off']], 80);
            uniform.set(value(v, 'uInteriorColor', [0, 0, 0, 0]), 84);
            uniform.set(value(v, 'uInteriorSubstance', [0, 0, 0, 0]), 88);
            uniform[92] = { off: 0, on: 1, opaque: 2 }[value<string>(v, 'dTransparentBackfaces', 'on') as 'on'];
            uniform[93] = value(v, 'uAlphaThickness', 0);
            uniform[94] = value(v, 'uDensity', 0.2);
            uniform[95] = value(v, 'dSolidInterior', false) ? 1 : 0;
            const grids = ['Substance', 'Overpaint', 'Transparency', 'Emissive', 'Color'].map((name, i) => {
                const type = value<string>(v, `d${name}Type`, 'groupInstance');
                const active = name === 'Color' ? type === 'volume' || type === 'volumeInstance' : value(v, `d${name}`, false) && type === 'volumeInstance';
                const data = value<unknown>(v, `t${name}Grid`, undefined);
                if (active && !(data instanceof WebGPUTextureData)) throw new Error(`Native spatial ${name} requires WebGPU-owned grid data.`);
                const native = active ? (data as WebGPUTextureData).native : undefined;
                if (native && native.device !== device) throw new Error('Native grid texture belongs to another GPU device.');
                const grid = native ? { width: native.texture.width, height: native.texture.height, depth: 1, array: new Uint8Array(0) } : active ? (data as WebGPUTextureData).data : { width: 1, height: 1, depth: 1, array: new Uint8Array(4) };
                if (grid.depth !== 1 || !(grid.array instanceof Uint8Array)) throw new Error(`Spatial ${name} requires a packed RGBA byte grid.`);
                uniform.set(value(v, `u${name}GridDim`, [1, 1, 1]), 96 + i * 8);
                uniform[99 + i * 8] = active ? (name === 'Color' ? (type === 'volume' ? 1 : 2) : value(v, `u${name}Strength`, 1)) : -1;
                uniform.set(value(v, `u${name}GridTransform`, [0, 0, 0, 1]), 100 + i * 8);
                return { ...grid, native: native?.texture, version: active ? this.valueVersions(v, [`t${name}Grid`]) : 'inactive' };
            });
            uniform[136] = value(v, 'dUsePalette', false) ? 1 : 0;
            uniform[137] = value(v, 'tPalette', { filter: 'nearest' }).filter === 'linear' ? 1 : 0;
            uniform[138] = value(v, 'uOverpaintStrength', 1);
            uniform[139] = value(v, 'dOverpaint', false) && value<string>(v, 'dOverpaintType', '') !== 'volumeInstance' ? 1 : 0;
            uniform.set(value(v, 'uLod', [0, 0, 0, 0]), 140);
            buffers.push(updateBuffer(4, uniform, GPUBufferUsage.UNIFORM));
            const image = geometry.kind === 3 ? value<{ width: number, height: number, array: Uint8Array } | undefined>(v, 'tFont', undefined) : geometry.kind === 4 ? value<{ width: number, height: number, array: Uint8Array } | undefined>(v, 'tImageTex', undefined) : undefined;
            const width = Math.max(1, image?.width || 1), height = Math.max(1, image?.height || 1);
            const single = geometry.kind === 3;
            texture = uploadTexture(0, geometry.kind === 3 || geometry.kind === 4 ? this.valueVersions(v, [single ? 'tFont' : 'tImageTex']) + `/${single}` : 'white',
                () => ({ width, height, array: image?.array || new Uint8Array([255, 255, 255, 255]) }), single ? 'r8unorm' : 'rgba8unorm');
            const group = value(v, 'tGroupTex', { array: new Uint8Array(4), width: 1, height: 1 });
            const scalar = value(v, 'tValueTex', { array: new Float32Array(1), width: 1, height: 1 });
            const palette = value(v, 'tPalette', { array: new Uint8Array([255, 255, 255]), width: 1, height: 1 });
            for (const [i, name] of ['tGroupTex', 'tValueTex', 'tPalette'].entries()) {
                const t = uploadTexture(i + 1, this.valueVersions(v, [name]), () => {
                    if (i === 0) return group;
                    if (i === 1) return scalar;
                    const rgba = new Uint8Array(palette.width * palette.height * 4);
                    for (let j = 0; j < rgba.length / 4; j++) { rgba.set(palette.array.subarray(j * 3, j * 3 + 3), j * 4); rgba[j * 4 + 3] = 255; }
                    return { ...palette, array: rgba };
                }, i === 1 ? 'r32float' : 'rgba8unorm');
                imageTextures.push(t); ownedImageTextures.push(t);
            }
            for (const [i, grid] of grids.entries()) {
                if (grid.native) { textureVersions[i + 4] = grid.version; imageTextures.push(grid.native); continue; }
                const t = uploadTexture(i + 4, grid.version, () => grid, 'rgba8unorm');
                imageTextures.push(t); ownedImageTextures.push(t);
            }
            buffers.push(updateBuffer(5, object.type === 'image' ? createImageStyles(v) : new Float32Array(8), GPUBufferUsage.STORAGE));
            buffers.push(updateBuffer(6, createClipData(v), GPUBufferUsage.STORAGE));
            buffers.push(updateBuffer(7, value(v, 'tClipping', { array: new Uint8Array(1) }).array, GPUBufferUsage.STORAGE));
            buffers.push(updateBuffer(8, value(v, 'tWiggle', { array: new Uint8Array(1) }).array, GPUBufferUsage.STORAGE));
            buffers.push(updateBuffer(9, createWebGPUOverpaintThemes(object, geometry), GPUBufferUsage.STORAGE));
            const bindGroup = device.createBindGroup({ layout: this.objectLayout, entries: [
                { binding: 0, resource: { buffer: buffers[4] } },
                ...[1, 2, 3].map((binding, index) => ({ binding, resource: { buffer: buffers[index] } })),
                { binding: 4, resource: texture.createView() },
                { binding: 5, resource: device.createSampler({ minFilter: 'linear', magFilter: 'linear' }) },
                ...imageTextures.slice(0, 3).map((t, i) => ({ binding: i + 6, resource: t.createView() })),
                ...imageTextures.slice(3).map((t, i) => ({ binding: i + 13, resource: t.createView() })),
                { binding: 9, resource: { buffer: buffers[5] } },
                { binding: 10, resource: { buffer: buffers[6] } },
                { binding: 11, resource: { buffer: buffers[7] } },
                { binding: 12, resource: { buffer: buffers[8] } },
                { binding: 18, resource: { buffer: buffers[9] } },
            ] });
            let transparent = geometry.kind === 3 || geometry.kind === 4 || (geometry.kind === 2 && value<string>(v, 'dPointStyle', 'square') === 'fuzzy') || value<string>(v, 'dXrayShaded', 'off') !== 'off';
            const spatialTransparency = value(v, 'dTransparency', false) && value<string>(v, 'dTransparencyType', '') === 'volumeInstance';
            transparent ||= spatialTransparency;
            let hasOpaque = spatialTransparency || value<string>(v, 'dTransparentBackfaces', 'off') === 'opaque';
            let themeHasTransparency = previous?.themeHasTransparency ?? false, themeHasOpaque = previous?.themeHasOpaque ?? false;
            if (themes) {
                themeHasTransparency = false; themeHasOpaque = false;
                for (let i = 3; i < themes.length; i += WebGPUThemeStride) { if (themes[i] < 1) themeHasTransparency = true; else themeHasOpaque = true; }
            }
            transparent ||= themeHasTransparency;
            hasOpaque ||= themeHasOpaque && value<string>(v, 'dXrayShaded', 'off') === 'off';
            const item = { object, buffers, borrowedGeometry: geometry.native?.buffer, texture, imageTextures, ownedImageTextures, bindGroup, count: geometry.indices.length, instances, version: this.version(object), transparent, hasOpaque, surface: geometry.kind === 0 || geometry.kind === 5 || geometry.kind === 6, geometry, geometryVersion, transformVersion, textureVersions, themeVersion, themeHasTransparency, themeHasOpaque };
            // Commit mutations only after validation and bind-group construction;
            // rejected updates leave the previous item's live buffers intact.
            for (const write of pendingWrites) {
                device.queue.writeBuffer(write.buffer, 0, write.data);
                this.bufferResourceStats.updates++; this.bufferResourceStats.uploadedBytes += write.data.byteLength;
            }
            if (!reuseThemes) this.themeResourceStats.uploads++;
            return item;
        } catch (error) {
            for (const buffer of allocatedBuffers) buffer.destroy();
            for (const t of allocatedTextures) t.destroy();
            throw error;
        }
    }

    private destroyItem(item: RenderItem, retained: readonly GPUBuffer[] = [], retainedTextures: readonly GPUTexture[] = []) {
        for (const buffer of item.buffers) if (buffer !== item.borrowedGeometry && !retained.includes(buffer)) buffer.destroy();
        if (!retainedTextures.includes(item.texture)) item.texture.destroy();
        for (const t of item.ownedImageTextures) if (!retainedTextures.includes(t)) t.destroy();
    }

    private resize(width: number, height: number) {
        if (width === this.width && height === this.height) return;
        if (width > this.context.device.limits.maxTextureDimension2D || height > this.context.device.limits.maxTextureDimension2D) throw new Error('WebGPU canvas exceeds the device texture dimension limit.');
        this.color?.destroy(); this.emissive?.destroy(); this.transparentColor?.destroy(); this.outlineDepthAlpha?.destroy(); this.outlineDepth?.destroy(); this.depth?.destroy(); this.colorPickDepth?.destroy(); this.pickTexture?.destroy(); this.opaqueDepth?.destroy(); this.selectionTexture?.destroy(); this.selectionDepth?.destroy(); this.selectionOpaqueDepth?.destroy();
        const size = { width, height };
        const usage = GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC | GPUTextureUsage.TEXTURE_BINDING;
        this.color = this.context.device.createTexture({ size, format: this.context.format, usage });
        this.emissive = this.context.device.createTexture({ size, format: 'rgba16float', usage });
        this.transparentColor = this.context.device.createTexture({ size, format: this.context.format, usage });
        this.outlineDepthAlpha = undefined; this.outlineDepth = undefined;
        this.depth = this.context.device.createTexture({ size, format: 'depth32float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC | GPUTextureUsage.TEXTURE_BINDING });
        this.colorPickDepth = this.context.device.createTexture({ size, format: 'depth32float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_DST });
        this.opaqueDepth = this.context.device.createTexture({ size, format: 'depth32float', usage: GPUTextureUsage.COPY_DST | GPUTextureUsage.TEXTURE_BINDING });
        this.pickTexture = this.context.device.createTexture({ size, format: 'rgba32uint', usage });
        this.selectionTexture = this.context.device.createTexture({ size, format: 'rgba32uint', usage: usage | GPUTextureUsage.COPY_DST });
        this.selectionDepth = this.context.device.createTexture({ size, format: 'depth32float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC });
        this.selectionOpaqueDepth = this.context.device.createTexture({ size, format: 'depth32float', usage: GPUTextureUsage.COPY_DST | GPUTextureUsage.TEXTURE_BINDING });
        this.width = width; this.height = height;
    }

    render(objects: readonly GraphicsRenderObject[], camera: Camera, props: RendererProps, transparentBackground = false, pixelRatio = 1, output?: { width: number, height: number, present: boolean }, postprocessing?: PostprocessingProps, marking?: MarkingProps, multiSample?: MultiSampleProps, changed = true, forceOn = false, illumination?: IlluminationProps) {
        this.stereoActive = false;
        this.outputPicking = this.selectionTexture;
        if (illumination?.enabled) {
            this.multiSample.reset();
            const key = JSON.stringify([illumination, multiSample, camera.viewport, output?.width ?? this.context.canvas.width, output?.height ?? this.context.canvas.height]);
            if (changed || key !== this.illuminationKey) this.illuminationIteration = 0;
            this.illuminationKey = key; this.illuminationMax = Math.pow(2, Math.max(0, Math.min(16, Math.round(illumination.maxIterations))));
            if (!this.illuminationNeedsFrame) return;
            // Tracing supplies opaque indirect lighting and directional shadows.
            // Keep transparent SSAO, but do not shade the opaque scene twice.
            const post = postprocessing ? { ...postprocessing, shadow: { name: 'off' as const, params: {} },
                outline: illumination.ignoreOutline ? { name: 'off' as const, params: {} } : postprocessing.outline } : undefined;
            if (multiSample?.mode === 'on') {
                if (this.disposed) throw new Error('WebGPU renderer has been disposed.');
                const width = output?.width ?? this.context.canvas.width, height = output?.height ?? this.context.canvas.height;
                if (!width || !height) return;
                this.resize(width, height);
                this.outputColor = this.multiSample.renderIllumination(camera, multiSample, this.illuminationIteration, this.illuminationMax, width, height, this.selectionTexture!, refresh => {
                    this.renderSingle(objects, camera, props, transparentBackground, pixelRatio, { width, height, present: false }, post, marking, undefined,
                        { props: illumination, iteration: this.illuminationIteration, refreshInput: refresh, fullSampling: true });
                    return this.outputColor!;
                });
                if (output?.present !== false && this.context.context) {
                    const encoder = this.context.device.createCommandEncoder();
                    encoder.copyTextureToTexture({ texture: this.outputColor }, { texture: this.context.context.getCurrentTexture() }, [width, height]); this.context.device.queue.submit([encoder.finish()]);
                }
            } else {
                this.renderSingle(objects, camera, props, transparentBackground, pixelRatio, output, post, marking, undefined, { props: illumination, iteration: this.illuminationIteration });
            }
            this.illuminationIteration++; return;
        }
        this.illuminationMax = 0; this.illuminationKey = '';
        if (!multiSample || multiSample.mode === 'off') { this.multiSample.reset(); this.renderSingle(objects, camera, props, transparentBackground, pixelRatio, output, postprocessing, marking); return; }
        if (this.disposed) throw new Error('WebGPU renderer has been disposed.');
        const width = output?.width ?? this.context.canvas.width, height = output?.height ?? this.context.canvas.height;
        if (!width || !height) return;
        this.resize(width, height);
        this.outputColor = this.multiSample.render(camera, multiSample, changed, forceOn, width, height, this.selectionTexture!, offset => {
            this.renderSingle(objects, camera, props, transparentBackground, pixelRatio, { width, height, present: false }, postprocessing, marking, offset);
            return this.outputColor!;
        });
        if (output?.present !== false && this.context.context) {
            const encoder = this.context.device.createCommandEncoder(); encoder.copyTextureToTexture({ texture: this.outputColor }, { texture: this.context.context.getCurrentTexture() }, [width, height]); this.context.device.queue.submit([encoder.finish()]);
        }
    }

    /** Independent eye accumulation, with canonical eye IDs/depth composed into one frame. */
    renderStereo(objects: readonly GraphicsRenderObject[], camera: Camera, stereoProps: StereoCameraProps, props: RendererProps, transparentBackground = false, pixelRatio = 1, output?: { width: number, height: number, present: boolean }, postprocessing?: PostprocessingProps, marking?: MarkingProps, multiSample?: MultiSampleProps, changed = true, forceOn = false) {
        if (camera.viewport.width < 2) { this.render(objects, camera, props, transparentBackground, pixelRatio, output, postprocessing, marking, multiSample, changed, forceOn); return; }
        if (this.disposed) throw new Error('WebGPU renderer has been disposed.');
        const { device, canvas } = this.context;
        const width = output?.width ?? canvas.width, height = output?.height ?? canvas.height;
        if (!width || !height) return;
        this.stereo ??= { camera: new WebGPUStereoCamera(), samples: [new WebGPUMultiSample(this.context), new WebGPUMultiSample(this.context)] };
        const stereo = this.stereo;
        const key = JSON.stringify([stereoProps, camera.viewport]);
        changed ||= stereo.key !== key; stereo.key = key;
        stereo.camera.update(camera, stereoProps);
        if (!stereo.color || stereo.color.width !== width || stereo.color.height !== height) {
            stereo.color?.destroy(); stereo.picking?.destroy();
            stereo.color = device.createTexture({ label: 'molstar-stereo-color', size: [width, height], format: this.context.format, usage: GPUTextureUsage.COPY_SRC | GPUTextureUsage.COPY_DST });
            stereo.picking = device.createTexture({ label: 'molstar-stereo-picking', size: [width, height], format: 'rgba32uint', usage: GPUTextureUsage.COPY_SRC | GPUTextureUsage.COPY_DST });
        }
        const total = { drawCount: 0, instanceCount: 0, triangleCount: 0 };
        try {
            for (const [i, eye] of [stereo.camera.left, stereo.camera.right].entries()) {
                this.multiSample = stereo.samples[i];
                // Shared postprocessing targets cannot reuse the other eye's AO history.
                const sampling = multiSample && { ...multiSample, reuseOcclusion: false };
                this.render(objects, eye, props, transparentBackground, pixelRatio, { width, height, present: false }, postprocessing, marking, sampling, changed, forceOn);
                for (const key of ['drawCount', 'instanceCount', 'triangleCount'] as const) total[key] += this.stats[key];
                const encoder = device.createCommandEncoder({ label: 'molstar-stereo-compose' });
                const origin = i === 0 ? { x: 0, y: 0 } : { x: eye.viewport.x, y: height - eye.viewport.y - eye.viewport.height };
                const size = i === 0 ? [width, height] : [eye.viewport.width, eye.viewport.height];
                encoder.copyTextureToTexture({ texture: this.outputColor!, origin }, { texture: stereo.color, origin }, size);
                encoder.copyTextureToTexture({ texture: this.selectionTexture!, origin }, { texture: stereo.picking!, origin }, size);
                device.queue.submit([encoder.finish()]);
            }
        } finally { this.multiSample = this.monoMultiSample; }
        this.stereoActive = true; this.outputColor = stereo.color; this.outputPicking = stereo.picking;
        Object.assign(this.stats, total);
        if (output?.present !== false && this.context.context) {
            const encoder = device.createCommandEncoder();
            encoder.copyTextureToTexture({ texture: stereo.color }, { texture: this.context.context.getCurrentTexture() }, [width, height]);
            device.queue.submit([encoder.finish()]);
        }
    }

    /** Prepare native illumination G-buffers without fog or environment composition. */
    renderTracingInput(objects: readonly GraphicsRenderObject[], camera: Camera, props: RendererProps, width = this.context.canvas.width, height = this.context.canvas.height) {
        this.renderSingle(objects, camera, props, true, 1, { width, height, present: false });
        const encoder = this.context.device.createCommandEncoder({ label: 'molstar-tracing-input' });
        this.drawTracingInput(encoder, camera, width, height);
        this.context.device.queue.submit([encoder.finish()]);
        return this.tracingInput.textures;
    }

    private drawTracingInput(encoder: GPUCommandEncoder, camera: Camera, width: number, height: number, backDepth = true) {
        this.tracingInput.render(encoder, camera, width, height, (pass, back) => {
            pass.setBindGroup(0, this.cameraBindGroup);
            for (const item of this.items.values()) {
                if (this.helperIds.has(item.object.id)) continue;
                if (!item.object.state.visible || !item.count || !item.instances || (item.surface ? !item.hasOpaque : item.transparent)) continue;
                pass.setPipeline(this.tracingInput.pipeline(back, cullBackfaces(item.object.values), value(item.object.values, 'dFlipSided', false), item.geometry.kind === 5));
                pass.setBindGroup(1, item.bindGroup); pass.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(pass, item, camera);
            }
        }, backDepth);
    }

    private renderSingle(objects: readonly GraphicsRenderObject[], camera: Camera, props: RendererProps, transparentBackground = false, pixelRatio = 1, output?: { width: number, height: number, present: boolean }, postprocessing?: PostprocessingProps, marking?: MarkingProps, occlusionOffset?: readonly number[], illumination?: { props: IlluminationProps, iteration: number, refreshInput?: boolean, fullSampling?: boolean }) {
        this.sphereLodRanges.clear(); this.sphereAnimation = props.enableAnimation;
        if (this.disposed) throw new Error('WebGPU renderer has been disposed.');
        const { canvas, device } = this.context;
        const width = output?.width ?? canvas.width, height = output?.height ?? canvas.height;
        if (!width || !height) return;
        this.resize(width, height);
        this.outputPicking = this.selectionTexture;
        const flattened = objects.flatMap(object => isMergedRenderObject(object) ? [...object.members] : [object]);
        const helper = this.cameraHelper?.isEnabled ? this.cameraHelper : undefined;
        helper?.update(camera);
        const helperObjects = helper?.getRenderObjects() ?? [];
        const handleObjects = this.handleHelper?.getRenderObjects() ?? [];
        const pointer = this.pointerHelper?.isEnabled ? this.pointerHelper : undefined;
        pointer?.setCamera(camera);
        const pointerObjects = pointer?.getRenderObjects() ?? [];
        const debugObjects = this.debugRegistry?.getRenderObjects() ?? [];
        const allHelpers = [...helperObjects, ...handleObjects, ...pointerObjects, ...debugObjects];
        this.helperIds.clear(); for (const object of allHelpers) this.helperIds.add(object.id);
        const volumes = flattened.filter(object => object.type === 'direct-volume');
        this.volumes.sync(volumes);
        const ids = new Set([...flattened, ...allHelpers].map(object => object.id));
        for (const [id, item] of this.items) {
            if (!ids.has(id)) { this.destroyItem(item); this.items.delete(id); }
        }
        const visible: RenderItem[] = [];
        for (const object of flattened) {
            if (object.type === 'direct-volume') continue;
            if (!object.state.visible || !value(object.values, 'drawCount', 0) || !value(object.values, 'instanceCount', 0)) continue;
            let item = this.items.get(object.id);
            if (!item || item.version !== this.version(object)) {
                const next = this.createItem(object, item);
                if (item) this.destroyItem(item, next.buffers, [next.texture, ...next.imageTextures]);
                this.items.set(object.id, next); item = next;
            }
            visible.push(item);
        }
        visible.sort((a, b) => {
            if (a.transparent !== b.transparent) return a.transparent ? 1 : -1;
            if (!a.transparent) return a.object.id - b.object.id;
            return Vec3.distance(camera.state.position, b.object.values.boundingSphere.ref.value.center) - Vec3.distance(camera.state.position, a.object.values.boundingSphere.ref.value.center);
        });
        const uniforms = new Float32Array(64);
        uniforms.set(camera.projection, 0); uniforms.set(getWebGPUModelView(camera), 16);
        const background = Color.toRgbNormalized(props.backgroundColor);
        const environmentReady = this.background.update(postprocessing?.enabled ? postprocessing.background : { variant: { name: 'off', params: {} } });
        const sceneTransparent = transparentBackground || environmentReady;
        const alpha = sceneTransparent ? 0 : 1;
        uniforms.set([...background, alpha], 32);
        uniforms.set([camera.viewport.width, camera.viewport.height, pixelRatio, camera.scale], 36);
        uniforms.set([camera.viewport.x, this.height - camera.viewport.y - camera.viewport.height, camera.viewport.width, camera.viewport.height], 52);
        uniforms.set([...Color.toRgbNormalized(props.highlightColor), props.highlightStrength], 40);
        uniforms.set([...Color.toRgbNormalized(props.selectColor), props.selectStrength], 44);
        uniforms.set([camera.state.fog ? camera.fogNear : 0, camera.state.fog ? camera.fogFar : 0, props.ambientIntensity, props.light[0]?.intensity ?? 1], 48);
        uniforms.set([this.time, props.enableAnimation ? 1 : 0, marking && marking.ghostEdgeStrength < 1 ? 1 : 0, environmentReady ? 1 : 0], 56);
        uniforms.set([height, camera.viewOffset.enabled ? camera.viewOffset.offsetX * 4 : 0, camera.viewOffset.enabled ? camera.viewOffset.offsetY * 4 : 0, 0], 60);
        device.queue.writeBuffer(this.cameraBuffer, 0, uniforms);
        device.queue.writeBuffer(this.lightingBuffer, 0, createLightingData(camera, props, width, height));
        const encoder = device.createCommandEncoder({ label: 'molstar-native-frame' });
        let pass = encoder.beginRenderPass({
            colorAttachments: [
                { view: this.color!.createView(), clearValue: { r: background[0] * alpha, g: background[1] * alpha, b: background[2] * alpha, a: alpha }, loadOp: 'clear', storeOp: 'store' },
                { view: this.pickTexture!.createView(), clearValue: { r: 0, g: 0, b: 0, a: 0 }, loadOp: 'clear', storeOp: 'store' },
                { view: this.emissive!.createView(), clearValue: { r: 0, g: 0, b: 0, a: 0 }, loadOp: 'clear', storeOp: 'store' },
            ],
            depthStencilAttachment: { view: this.depth!.createView(), depthClearValue: 1, depthLoadOp: 'clear', depthStoreOp: 'store' },
        });
        const viewport = camera.viewport;
        pass.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
        pass.setBindGroup(0, this.cameraBindGroup);
        this.stats.drawCount = 0; this.stats.instanceCount = 0; this.stats.triangleCount = 0;
        const emitTransparent = !(postprocessing?.bloom.name === 'on' && postprocessing.bloom.params.mode === 'emissive' && !postprocessing.bloom.params.transparency);
        const draw = (item: RenderItem, transparent: boolean, phase: 'all' | 'opaque' | 'transparent') => {
            pass.setPipeline(this.pipeline(transparent, cullBackfaces(item.object.values), !transparent || emitTransparent, value(item.object.values, 'dFlipSided', false), phase, item.geometry.kind === 5)); pass.setBindGroup(1, item.bindGroup);
            pass.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(pass, item, camera, true);
        };
        // Surface opacity can vary within an object; opaque fragments must populate depth first.
        for (const item of visible) if ((item.surface && item.hasOpaque) || !item.transparent) draw(item, false, item.surface ? 'opaque' : 'all');
        if (illumination) {
            pass.end();
            if (illumination.iteration === 0 || illumination.refreshInput) this.drawTracingInput(encoder, camera, width, height, illumination.props.thicknessMode === 'auto');
            const input = this.tracingInput.textures;
            this.quality = this.tracingQuality.update(illumination.props, illumination.iteration);
            const traced = this.tracing.render(camera, props, illumination.props, input, illumination.iteration, this.quality.rendersPerFrame, this.quality.steps, this.quality.refineSteps, encoder);
            this.illuminationCompose.render(encoder, traced, input, camera, illumination.props, props, sceneTransparent, illumination.iteration, illumination.fullSampling, this.color);
            pass = encoder.beginRenderPass({ label: 'molstar-illumination-transparent',
                colorAttachments: [this.color!, this.pickTexture!, this.emissive!].map(texture => ({ view: texture.createView(), loadOp: 'load' as const, storeOp: 'store' as const })),
                depthStencilAttachment: { view: this.depth!.createView(), depthLoadOp: 'load', depthStoreOp: 'store' },
            });
            pass.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
            pass.setBindGroup(0, this.cameraBindGroup);
        }
        const weighted = this.transparency === 'wboit';
        const ordered = this.transparency === 'dpoit', oit = weighted || ordered;
        if (!oit) for (const item of visible) if (item.transparent) draw(item, true, item.surface ? 'transparent' : 'all');
        pass.end();
        if (oit) {
            if (volumes.length || ordered) encoder.copyTextureToTexture({ texture: this.depth! }, { texture: this.opaqueDepth! }, [this.width, this.height]);
            if (weighted) {
                const accum = this.weightedTransparency.begin(encoder, this.depth!, camera);
                accum.setBindGroup(0, this.cameraBindGroup);
                for (const item of visible) if (item.transparent) {
                    accum.setPipeline(this.weightedPipeline(cullBackfaces(item.object.values), value(item.object.values, 'dFlipSided', false), emitTransparent, item.geometry.kind === 5));
                    accum.setBindGroup(1, item.bindGroup); accum.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(accum, item, camera, true);
                }
                if (volumes.length) this.stats.drawCount += this.volumes.render(accum, volumes, camera, this.cameraBindGroup, this.opaqueDepth!, emitTransparent, false, false, false, true);
                accum.end();
                this.weightedTransparency.resolve(encoder, this.color!, this.emissive!, camera, this.transparentColor!, emitTransparent);
            } else {
                this.depthPeeling.render(encoder, this.opaqueDepth!, this.color!, this.emissive!, this.transparentColor!, camera, this.peelingIterations, emitTransparent, (pass, phase, bindings) => {
                    pass.setBindGroup(0, this.cameraBindGroup); pass.setBindGroup(2, this.depthPeeling.emptyGroup); pass.setBindGroup(3, bindings);
                    for (const item of visible) if (item.transparent) {
                        pass.setPipeline(this.peelingPipeline(phase, cullBackfaces(item.object.values), value(item.object.values, 'dFlipSided', false), item.geometry.kind === 5));
                        pass.setBindGroup(1, item.bindGroup); pass.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(pass, item, camera, true);
                    }
                    if (volumes.length) this.stats.drawCount += this.volumes.render(pass, volumes, camera, this.cameraBindGroup, this.opaqueDepth!, emitTransparent, false, false, false, false, false, phase);
                });
            }
            // Postprocessing needs visible color depth even for faint or unpickable geometry.
            encoder.copyTextureToTexture({ texture: this.depth! }, { texture: this.colorPickDepth! }, [this.width, this.height]);
            const ids = encoder.beginRenderPass({ label: 'molstar-weighted-color-depth',
                colorAttachments: [{ view: this.pickTexture!.createView(), loadOp: 'load', storeOp: 'store' }],
                depthStencilAttachment: { view: this.colorPickDepth!.createView(), depthLoadOp: 'load', depthStoreOp: 'store' },
            });
            ids.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1); ids.setBindGroup(0, this.cameraBindGroup);
            for (const item of visible) if (item.transparent) {
                ids.setPipeline(this.colorDepthPipeline(cullBackfaces(item.object.values), value(item.object.values, 'dFlipSided', false), item.geometry.kind === 5));
                ids.setBindGroup(1, item.bindGroup); ids.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(ids, item, camera);
            }
            if (volumes.length) this.volumes.render(ids, volumes, camera, this.cameraBindGroup, this.opaqueDepth!, true, false, false, false, false, true);
            ids.end();
        } else if (volumes.length) {
            encoder.copyTextureToTexture({ texture: this.depth! }, { texture: this.opaqueDepth! }, { width: this.width, height: this.height });
            const volumePass = encoder.beginRenderPass({
                colorAttachments: [
                    { view: this.color!.createView(), loadOp: 'load', storeOp: 'store' },
                    { view: this.pickTexture!.createView(), loadOp: 'load', storeOp: 'store' },
                    { view: this.emissive!.createView(), loadOp: 'load', storeOp: 'store' },
                ],
                depthStencilAttachment: { view: this.depth!.createView(), depthReadOnly: true },
            });
            volumePass.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
            this.stats.drawCount += this.volumes.render(volumePass, volumes, camera, this.cameraBindGroup, this.opaqueDepth!, emitTransparent);
            volumePass.end();
        }
        let transparencyMin = 1;
        for (const item of visible) {
            const values = item.object.values;
            const alpha = Math.max(0, Math.min(1, value(values, 'alpha', 1) * item.object.state.alphaFactor));
            if (alpha < 1) transparencyMin = Math.min(transparencyMin, 1 - alpha);
            if (value<string>(values, 'dXrayShaded', 'off') !== 'off' || value<string>(values, 'dPointStyle', '') === 'fuzzy' || item.object.type === 'text' || item.object.type === 'image') transparencyMin = Math.min(transparencyMin, 0.5);
            const min = value(values, 'transparencyMin', 0);
            if (min > 0) transparencyMin = Math.min(transparencyMin, min);
        }
        const transparentAo = !!(postprocessing?.enabled && postprocessing.occlusion.name === 'on' && transparencyMin < postprocessing.occlusion.params.transparentThreshold);
        // Separate color pass respects the portable 32-byte MRT limit.
        if (!oit && postprocessing?.enabled && (postprocessing.outline.name === 'on' || postprocessing.occlusion.name === 'on' || postprocessing.shadow.name === 'on')) {
            const transparentPass = encoder.beginRenderPass({ label: 'molstar-transparent-outline-color',
                colorAttachments: [{ view: this.transparentColor!.createView(), loadOp: 'clear', storeOp: 'store' }],
                depthStencilAttachment: { view: this.depth!.createView(), depthReadOnly: true },
            });
            transparentPass.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
            transparentPass.setBindGroup(0, this.cameraBindGroup);
            for (const item of visible) if (item.transparent) {
                transparentPass.setPipeline(this.transparentColorPipeline(cullBackfaces(item.object.values), value(item.object.values, 'dFlipSided', false), item.surface, item.geometry.kind === 5));
                transparentPass.setBindGroup(1, item.bindGroup); transparentPass.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(transparentPass, item, camera);
            }
            if (volumes.length) this.volumes.render(transparentPass, volumes, camera, this.cameraBindGroup, this.opaqueDepth!, false, false, false, true);
            transparentPass.end();
        }
        if (transparentAo || (postprocessing?.enabled && postprocessing.outline.name === 'on' && postprocessing.outline.params.includeTransparent)) {
            if (!this.outlineDepthAlpha) {
                const size = [this.width, this.height];
                this.outlineDepthAlpha = device.createTexture({ size, format: 'rgba32float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING });
                this.outlineDepth = device.createTexture({ size, format: 'depth32float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_DST });
            }
            encoder.copyTextureToTexture({ texture: this.depth! }, { texture: this.outlineDepth! }, [this.width, this.height]);
            const outlinePass = encoder.beginRenderPass({ label: 'molstar-transparent-outline-depth',
                colorAttachments: [{ view: this.outlineDepthAlpha.createView(), loadOp: 'clear', storeOp: 'store', clearValue: { r: 1, g: 0, b: 0, a: 0 } }],
                depthStencilAttachment: { view: this.outlineDepth!.createView(), depthLoadOp: 'load', depthStoreOp: 'store' },
            });
            outlinePass.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
            outlinePass.setBindGroup(0, this.cameraBindGroup);
            for (const item of visible) if (item.transparent) {
                outlinePass.setPipeline(this.outlineDepthPipeline(cullBackfaces(item.object.values), value(item.object.values, 'dFlipSided', false), item.surface, item.geometry.kind === 5));
                outlinePass.setBindGroup(1, item.bindGroup); outlinePass.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(outlinePass, item, camera);
            }
            if (volumes.length) this.volumes.render(outlinePass, volumes, camera, this.cameraBindGroup, this.opaqueDepth!, false, false, true);
            outlinePass.end();
        }
        // Selection has its own depth and IDs: rejected faint fragments cannot erase geometry behind them.
        const pickingPass = encoder.beginRenderPass({ label: 'molstar-picking',
            colorAttachments: [{ view: this.selectionTexture!.createView(), loadOp: 'clear', storeOp: 'store' }],
            depthStencilAttachment: { view: this.selectionDepth!.createView(), depthClearValue: 1, depthLoadOp: 'clear', depthStoreOp: 'store' },
        });
        pickingPass.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
        pickingPass.setBindGroup(0, this.cameraBindGroup);
        for (const item of visible) {
            if (!item.object.state.pickable || item.object.state.colorOnly) continue;
            pickingPass.setPipeline(this.pickingPipeline(cullBackfaces(item.object.values), value(item.object.values, 'dFlipSided', false), item.geometry.kind === 5));
            pickingPass.setBindGroup(1, item.bindGroup); pickingPass.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(pickingPass, item, camera);
        }
        pickingPass.end();
        if (volumes.length) {
            encoder.copyTextureToTexture({ texture: this.selectionDepth! }, { texture: this.selectionOpaqueDepth! }, { width: this.width, height: this.height });
            const pickingVolumes = encoder.beginRenderPass({ label: 'molstar-volume-picking',
                colorAttachments: [{ view: this.selectionTexture!.createView(), loadOp: 'load', storeOp: 'store' }],
                depthStencilAttachment: { view: this.selectionDepth!.createView(), depthLoadOp: 'load', depthStoreOp: 'store' },
            });
            pickingVolumes.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
            this.volumes.render(pickingVolumes, volumes, camera, this.cameraBindGroup, this.selectionOpaqueDepth!, false, true);
            pickingVolumes.end();
        }
        let sceneCenter = camera.state.target;
        if (postprocessing?.enabled && postprocessing.dof.name === 'on' && postprocessing.dof.params.center === 'scene-center') {
            const spheres = flattened.filter(o => o.state.visible).map(o => o.values.boundingSphere.ref.value);
            if (spheres.length) {
                const boundary = new BoundaryHelper('98');
                for (const sphere of spheres) boundary.includeSphere(sphere);
                boundary.finishedIncludeStep();
                for (const sphere of spheres) boundary.radiusSphere(sphere);
                sceneCenter = boundary.getSphere().center;
            }
        }
        const hasEmission = (o: GraphicsRenderObject) => o.state.visible && (value(o.values, 'uEmissive', 0) > 0 || value(o.values, 'dEmissive', false));
        const emits = visible.some(item => (item.hasOpaque || !item.transparent || emitTransparent) && hasEmission(item.object)) || (emitTransparent && volumes.some(hasEmission));
        let applyMarking: ((source: GPUTexture) => GPUTexture) | undefined;
        if (marking?.enabled && visible.some(item => value(item.object.values, 'markerAverage', 0) > 0 || value(item.object.values, 'uMarker', 0) > 0)) {
            this.marking.setSize(this.width, this.height);
            const drawMarkers = (mask: boolean) => {
                const pass = encoder.beginRenderPass({ label: mask ? 'molstar-marked-mask' : 'molstar-unmarked-depth',
                    colorAttachments: mask ? [{ view: this.marking.mask!.createView(), clearValue: { r: 1, g: 1, b: 1, a: 1 }, loadOp: 'clear', storeOp: 'store' }] : [],
                    depthStencilAttachment: { view: (mask ? this.marking.maskDepth! : this.marking.depth!).createView(), depthClearValue: 1, depthLoadOp: 'clear', depthStoreOp: 'store' },
                });
                pass.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
                pass.setBindGroup(0, this.cameraBindGroup);
                for (const item of visible) {
                    if (mask && value(item.object.values, 'markerAverage', 0) <= 0 && value(item.object.values, 'uMarker', 0) <= 0) continue;
                    pass.setPipeline(this.markingPipeline(mask, cullBackfaces(item.object.values), value(item.object.values, 'dFlipSided', false), item.geometry.kind === 5));
                    if (mask) pass.setBindGroup(2, this.marking.geometryBindings());
                    pass.setBindGroup(1, item.bindGroup); pass.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(pass, item, camera);
                }
                pass.end();
            };
            if (marking.ghostEdgeStrength < 1) drawMarkers(false);
            drawMarkers(true);
            applyMarking = source => this.marking.render(encoder, source, camera, marking, pixelRatio);
        }
        const prepareHelpers = (objects: GraphicsRenderObject[]) => objects.filter(object => object.state.visible && value(object.values, 'drawCount', 0) > 0 && value(object.values, 'instanceCount', 0) > 0).map(object => {
            let item = this.items.get(object.id);
            if (!item || item.version !== this.version(object)) {
                const next = this.createItem(object, item);
                if (item) this.destroyItem(item, next.buffers, [next.texture, ...next.imageTextures]);
                this.items.set(object.id, next); item = next;
            }
            return item;
        });
        const drawWorldHelpers = (source: GPUTexture, objects: GraphicsRenderObject[], bindings: GPUBindGroup, pickable: boolean, name: string) => {
            const items = prepareHelpers(objects);
            if (!items.length) return source;
            const pass = encoder.beginRenderPass({ label: `molstar-${name}-helper-color`,
                colorAttachments: [source, this.pickTexture!, this.emissive!].map(texture => ({ view: texture.createView(), loadOp: 'load' as const, storeOp: 'store' as const })),
                depthStencilAttachment: { view: this.depth!.createView(), depthLoadOp: 'load', depthStoreOp: 'store' },
            });
            pass.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1); pass.setBindGroup(0, bindings);
            for (const item of items) {
                pass.setPipeline(this.pipeline(item.transparent, cullBackfaces(item.object.values), true, false, 'all', item.geometry.kind === 5));
                pass.setBindGroup(1, item.bindGroup); pass.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(pass, item, camera, true);
            }
            pass.end();
            if (pickable) {
                const pick = encoder.beginRenderPass({ label: `molstar-${name}-helper-picking`,
                    colorAttachments: [{ view: this.selectionTexture!.createView(), loadOp: 'load', storeOp: 'store' }],
                    depthStencilAttachment: { view: this.selectionDepth!.createView(), depthLoadOp: 'load', depthStoreOp: 'store' },
                });
                pick.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1); pick.setBindGroup(0, bindings);
                for (const item of items) {
                    if (!item.object.state.pickable || item.object.state.colorOnly) continue;
                    pick.setPipeline(this.pickingPipeline(cullBackfaces(item.object.values), value(item.object.values, 'dFlipSided', false), item.geometry.kind === 5));
                    pick.setBindGroup(1, item.bindGroup); pick.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(pick, item, camera);
                }
                pick.end();
            }
            return source;
        };
        if (pointer) {
            const pointerUniforms = new Float32Array(uniforms);
            pointerUniforms.set(pointer.camera.projection, 0); pointerUniforms.set(getWebGPUModelView(pointer.camera), 16); pointerUniforms[39] = 1;
            device.queue.writeBuffer(this.pointerCameraBuffer, 0, pointerUniforms);
        }
        const drawHelper = (source: GPUTexture) => {
            if (!helper || !helperObjects.length) return source;
            if (!this.helperDepth || this.helperDepth.width !== width || this.helperDepth.height !== height) {
                this.helperDepth?.destroy();
                this.helperDepth = device.createTexture({ label: 'molstar-camera-helper-depth', size: [width, height], format: 'depth32float', usage: GPUTextureUsage.RENDER_ATTACHMENT });
            }
            const helperUniforms = new Float32Array(uniforms);
            helperUniforms.set(helper.camera.projection, 0);
            helperUniforms.set(Mat4.mul(Mat4(), helper.camera.view, helper.scene.view), 16);
            helperUniforms[39] = 1; // Helpers are specified in display pixels, independent of model scale.
            helperUniforms[48] = 0; helperUniforms[49] = 0; // Orientation axes never inherit molecular fog.
            helperUniforms.set([0, 0, 0, 0], 56);
            device.queue.writeBuffer(this.helperCameraBuffer, 0, helperUniforms);
            const helperItems = prepareHelpers(helperObjects);
            const overlay = encoder.beginRenderPass({ label: 'molstar-camera-helper-color',
                colorAttachments: [source, this.pickTexture!, this.emissive!].map(texture => ({ view: texture.createView(), loadOp: 'load' as const, storeOp: 'store' as const })),
                depthStencilAttachment: { view: this.helperDepth.createView(), depthLoadOp: 'clear', depthClearValue: 1, depthStoreOp: 'store' },
            });
            overlay.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
            overlay.setBindGroup(0, this.helperCameraBindGroup);
            for (const item of helperItems) {
                overlay.setPipeline(this.pipeline(item.transparent, cullBackfaces(item.object.values), true, false, 'all', item.geometry.kind === 5));
                overlay.setBindGroup(1, item.bindGroup); overlay.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(overlay, item, camera, true);
            }
            overlay.end();
            const picking = encoder.beginRenderPass({ label: 'molstar-camera-helper-picking',
                colorAttachments: [{ view: this.selectionTexture!.createView(), loadOp: 'load', storeOp: 'store' }],
                depthStencilAttachment: { view: this.helperDepth.createView(), depthLoadOp: 'clear', depthClearValue: 1, depthStoreOp: 'store' },
            });
            picking.setViewport(viewport.x, this.height - viewport.y - viewport.height, viewport.width, viewport.height, 0, 1);
            picking.setBindGroup(0, this.helperCameraBindGroup);
            for (const item of helperItems) {
                picking.setPipeline(this.pickingPipeline(cullBackfaces(item.object.values), false, item.geometry.kind === 5));
                picking.setBindGroup(1, item.bindGroup); picking.setIndexBuffer(item.buffers[3], 'uint32'); this.drawItem(picking, item, camera);
            }
            picking.end();
            return source;
        };
        const finishOverlays = allHelpers.length ? (source: GPUTexture) => {
            source = applyMarking?.(source) ?? source;
            source = drawWorldHelpers(source, debugObjects, this.cameraBindGroup, false, 'debug');
            source = drawWorldHelpers(source, pointerObjects, this.pointerCameraBindGroup, false, 'pointer');
            source = drawWorldHelpers(source, handleObjects, this.cameraBindGroup, true, 'handle');
            return drawHelper(source);
        } : applyMarking;
        this.outputColor = this.postprocessing.render(encoder, this.color!, this.depth!, this.pickTexture!, camera, postprocessing, pixelRatio, sceneCenter, this.emissive!, emits, props, sceneTransparent, this.outlineDepthAlpha, this.transparentColor, transparentAo, finishOverlays, environmentReady ? source => this.background.render(encoder, source, camera, postprocessing!.background, transparentBackground, props.backgroundColor) : undefined, occlusionOffset, !illumination);
        if (output?.present !== false && this.context.context) encoder.copyTextureToTexture({ texture: this.outputColor }, { texture: this.context.context.getCurrentTexture() }, { width: this.width, height: this.height });
        device.queue.submit([encoder.finish()]);
    }

    async pick(x: number, y: number): Promise<WebGPUPick | undefined> {
        if (!this.outputPicking || x < 0 || y < 0 || x >= this.width || y >= this.height) return undefined;
        const pixels = await this.context.readTexture(this.outputPicking, Math.floor(x), Math.floor(y), 1, 1, 16);
        const ids = new Uint32Array(pixels.buffer);
        if (!ids[0]) return undefined;
        const object = this.items.get(ids[0] - 1)?.object ?? this.volumes.getObject(ids[0] - 1);
        if (!object?.state.visible || !object.state.pickable) return undefined;
        return { id: { objectId: ids[0] - 1, instanceId: ids[1], groupId: ids[2] }, depth: new Float32Array(pixels.buffer)[3] };
    }

    async readPixels(): Promise<{ width: number, height: number, array: Uint8Array }> {
        if (!this.outputColor) throw new Error('Render a WebGPU frame before reading pixels.');
        const width = this.width, height = this.height;
        const array = await this.context.readTexture(this.outputColor, 0, 0, width, height);
        if (this.context.format.startsWith('bgra')) {
            for (let i = 0; i < array.length; i += 4) { const r = array[i]; array[i] = array[i + 2]; array[i + 2] = r; }
        }
        return { width, height, array };
    }

    /** Canonical packed object/instance/group IDs and float32 depth bits. */
    async readPickingPixels(): Promise<{ width: number, height: number, array: Uint32Array }> {
        if (!this.outputPicking) throw new Error('Render a WebGPU frame before reading picking pixels.');
        const width = this.width, height = this.height;
        const bytes = await this.context.readTexture(this.outputPicking, 0, 0, width, height, 16);
        return { width, height, array: new Uint32Array(bytes.buffer, bytes.byteOffset, bytes.byteLength / 4) };
    }

    dispose() {
        this.sphereLodRanges.clear();
        if (this.disposed) return;
        this.disposed = true;
        for (const item of this.items.values()) this.destroyItem(item);
        this.items.clear(); this.pipelines.clear();
        this.volumes.dispose(); this.postprocessing.dispose(); this.marking.dispose(); this.background.dispose(); this.multiSample.dispose(); this.tracingInput.dispose(); this.tracing.dispose(); this.illuminationCompose.dispose();
        this.weightedTransparency.dispose();
        this.depthPeeling.dispose();
        if (this.stereo) { for (const sample of this.stereo.samples) sample.dispose(); this.stereo.color?.destroy(); this.stereo.picking?.destroy(); }
        this.cameraBuffer.destroy(); this.helperCameraBuffer.destroy(); this.pointerCameraBuffer.destroy(); this.helperDepth?.destroy(); this.lightingBuffer.destroy(); this.color?.destroy(); this.emissive?.destroy(); this.transparentColor?.destroy(); this.outlineDepthAlpha?.destroy(); this.outlineDepth?.destroy(); this.depth?.destroy(); this.colorPickDepth?.destroy(); this.pickTexture?.destroy(); this.opaqueDepth?.destroy(); this.selectionTexture?.destroy(); this.selectionDepth?.destroy(); this.selectionOpaqueDepth?.destroy();
    }
}

function cullBackfaces(values: RenderableValues) {
    return !value(values, 'uDoubleSided', true) && !value(values, 'dSolidInterior', false);
}
