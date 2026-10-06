/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { getLight, getTransformedLightDirection, RendererProps } from '../renderer';
import { WebGPUSmaa } from './smaa';
import { WebGPUBloom } from './bloom';
import { WebGPUDepthPyramid } from './depth-pyramid';
import { Camera } from '../../mol-canvas3d/camera';
import { getSsaoSamples } from '../../mol-canvas3d/passes/ssao';
import { PostprocessingProps } from '../../mol-canvas3d/passes/postprocessing';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Color } from '../../mol-util/color';
import { WebGPUContext } from './context';
import { GPUBufferUsage, GPUShaderStage, GPUTextureUsage } from './compat';

/** Native screen-space passes. Picking remains attached to the original geometry. */
export class WebGPUPostprocessing {
    private readonly bloom: WebGPUBloom;
    private readonly depthPyramid: WebGPUDepthPyramid;
    private readonly uniforms = new Map<string, GPUBuffer>();
    private samples: GPUBuffer;
    private levels: GPUBuffer;
    private lights: GPUBuffer;
    private lightsKey = '';
    private sampleCount = 0;
    private levelsKey = '';
    private outlineTarget?: GPUTexture;
    private transparentAoColor?: GPUTexture;
    private readonly aoTargets: GPUTexture[] = [];
    private readonly layout: GPUBindGroupLayout;
    private readonly sampler: GPUSampler;
    private readonly pipelines: Record<'outline' | 'fxaa' | 'ssao' | 'ssaoBlur' | 'ssaoCompose' | 'ssaoTransparentColor' | 'outlineEdges' | 'sharpen' | 'dof' | 'bloomCompose' | 'shadowCompose', GPURenderPipeline>;
    private readonly targets: GPUTexture[] = [];
    private width = 0;
    private height = 0;
    private aoWidth = 0;
    private aoHeight = 0;

    constructor(private readonly context: WebGPUContext, shader: GPUShaderModule, bloomModule: GPUShaderModule, private readonly smaa: WebGPUSmaa) {
        this.bloom = new WebGPUBloom(context, bloomModule);
        this.depthPyramid = new WebGPUDepthPyramid(context);
        const { device } = context;
        this.lights = context.createBuffer(new Float32Array(8), GPUBufferUsage.STORAGE);
        this.samples = context.createBuffer(new Float32Array(4), GPUBufferUsage.STORAGE);
        this.levels = context.createBuffer(new Float32Array(4), GPUBufferUsage.STORAGE);
        this.layout = device.createBindGroupLayout({ entries: [
            { binding: 0, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
            { binding: 1, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 2, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
            { binding: 3, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'depth' } },
            { binding: 4, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'uint' } },
            { binding: 5, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            ...[8, 9].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' as const } })),
            { binding: 11, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 12, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'unfilterable-float' } },
            { binding: 13, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
            { binding: 14, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'unfilterable-float' } },
            ...[6, 7, 10].map(binding => ({ binding, visibility: GPUShaderStage.FRAGMENT, buffer: { type: 'read-only-storage' as const } })),
        ] });
        const layout = device.createPipelineLayout({ bindGroupLayouts: [this.layout] });
        const pipeline = (entryPoint: string, format: GPUTextureFormat = context.format) => device.createRenderPipeline({ label: `molstar-${entryPoint}`, layout,
            vertex: { module: shader, entryPoint: 'vs' }, fragment: { module: shader, entryPoint, targets: [{ format }] }, primitive: { topology: 'triangle-list' } });
        this.pipelines = { outline: pipeline('outline'), fxaa: pipeline('fxaa'), ssao: pipeline('ssao', 'rgba16float'), ssaoBlur: pipeline('ssaoBlur', 'rgba16float'), ssaoCompose: pipeline('ssaoCompose'), ssaoTransparentColor: pipeline('ssaoTransparentColor'), outlineEdges: pipeline('outlineEdges', 'rgba16float'), sharpen: pipeline('sharpen'), dof: pipeline('dof'), bloomCompose: pipeline('bloomCompose'), shadowCompose: pipeline('shadowCompose') };
        this.sampler = device.createSampler({ minFilter: 'linear', magFilter: 'linear' });
    }

    render(encoder: GPUCommandEncoder, color: GPUTexture, depth: GPUTexture, picking: GPUTexture, camera: Camera, props: PostprocessingProps | undefined, pixelRatio: number, sceneCenter = camera.state.target, emissive = color, hasEmission = true, rendererProps?: RendererProps, transparentBackground = true, transparentDepth = emissive, transparentColor = emissive, includeTransparentAo = false, applyMarking?: (source: GPUTexture) => GPUTexture, applyBackground?: (source: GPUTexture) => GPUTexture, occlusionOffset?: readonly number[], includeOpaqueAo = true) {
        const finish = (source: GPUTexture) => { const composited = applyBackground?.(source) ?? source; return applyMarking?.(composited) ?? composited; };
        if (!props?.enabled) return finish(color);
        const stages: ('outline' | 'fxaa' | 'ssaoCompose' | 'sharpen' | 'dof' | 'bloomCompose' | 'shadowCompose' | 'smaa')[] = [];
        const aoEnabled = props.occlusion.name === 'on' && (includeOpaqueAo || includeTransparentAo);
        if (aoEnabled) stages.push('ssaoCompose');
        const shadowEnabled = props.shadow.name === 'on' && !!rendererProps;
        if (shadowEnabled) stages.push('shadowCompose');
        if (props.outline.name === 'on') stages.push('outline');
        const bloomEnabled = props.bloom.name === 'on' && props.bloom.params.strength > 0 && (props.bloom.params.mode === 'luminosity' || hasEmission);
        if (bloomEnabled) stages.push('bloomCompose');
        if (props.dof.name === 'on') stages.push('dof');
        if (props.antialiasing.name === 'fxaa') stages.push('fxaa');
        if (props.antialiasing.name === 'smaa') stages.push('smaa');
        if (props.sharpening.name === 'on') stages.push('sharpen');
        if (!stages.length) return finish(color);
        if (color.width !== this.width || color.height !== this.height) {
            for (const t of [...this.targets, ...this.aoTargets]) t.destroy(); this.targets.length = 0; this.aoTargets.length = 0; this.outlineTarget?.destroy(); this.outlineTarget = undefined; this.transparentAoColor?.destroy(); this.transparentAoColor = undefined;
            this.width = color.width; this.height = color.height;
            for (let i = 0; i < 2; i++) this.targets.push(this.context.device.createTexture({ size: { width: this.width, height: this.height }, format: this.context.format,
                usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_SRC }));
        }
        const aoScale = props.occlusion.name === 'on' ? Math.min(1, 1 / pixelRatio) * props.occlusion.params.resolutionScale : 1;
        const aoWidth = Math.max(1, Math.floor(this.width * aoScale)), aoHeight = Math.max(1, Math.floor(this.height * aoScale));
        const reuseAo = !!occlusionOffset && this.aoTargets.length > 0 && aoWidth === this.aoWidth && aoHeight === this.aoHeight;
        if (aoEnabled && (!this.aoTargets.length || aoWidth !== this.aoWidth || aoHeight !== this.aoHeight)) {
            for (const t of this.aoTargets) t.destroy(); this.aoTargets.length = 0;
            this.aoWidth = aoWidth; this.aoHeight = aoHeight;
            for (let i = 0; i < 2; i++) this.aoTargets.push(this.context.device.createTexture({ label: 'molstar-scaled-ssao', size: { width: aoWidth, height: aoHeight }, format: 'rgba16float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING }));
        }
        if (aoEnabled && includeTransparentAo && !this.transparentAoColor) this.transparentAoColor = this.context.device.createTexture({ size: { width: this.width, height: this.height }, format: this.context.format, usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING });
        if (props.outline.name === 'on' && !this.outlineTarget) this.outlineTarget = this.context.device.createTexture({ size: { width: this.width, height: this.height }, format: 'rgba16float', usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.TEXTURE_BINDING });
        const settings = new Float32Array(96);
        settings.set([aoWidth, aoHeight, reuseAo ? occlusionOffset![0] / this.width : 0, reuseAo ? occlusionOffset![1] / this.height : 0], 92);
        settings.set(rendererProps ? Color.toRgbNormalized(rendererProps.backgroundColor) : [0, 0, 0], 88);
        settings[87] = camera.state.fog ? camera.fogNear : 0;
        settings[91] = camera.state.fog ? camera.fogFar : 0;
        settings.set(Mat4.invert(Mat4(), camera.projection));
        settings.set([this.width, this.height, includeTransparentAo ? 1 : 0, camera.far], 16);
        settings[34] = transparentBackground ? 1 : 0;
        const v = camera.viewport;
        settings.set([v.x, this.height - v.y - v.height, v.width, v.height], 20);
        if (props.antialiasing.name === 'fxaa') {
            const p = props.antialiasing.params;
            settings.set([p.edgeThresholdMin, p.edgeThresholdMax, p.iterations, p.subpixelQuality], 24);
        }
        if (props.outline.name === 'on') {
            const p = props.outline.params;
            settings.set([...Color.toRgbNormalized(p.color), p.threshold], 28);
            settings.set([Math.max(1, Math.round(p.scale * pixelRatio)), p.includeTransparent ? 1 : 0, transparentBackground ? 1 : 0, 0], 32);
        }
        settings[35] = applyBackground ? 1 : 0;
        settings.set(camera.projection, 36);
        if (props.occlusion.name === 'on') {
            const p = props.occlusion.params;
            if (this.sampleCount !== p.samples) {
                const packed = new Float32Array(p.samples * 4); const samples = getSsaoSamples(p.samples);
                for (let i = 0; i < p.samples; i++) packed.set(samples.slice(i * 3, i * 3 + 3), i * 4);
                this.samples.destroy(); this.samples = this.context.createBuffer(packed, GPUBufferUsage.STORAGE); this.sampleCount = p.samples;
            }
            const multi = p.multiScale;
            const levels = multi.name === 'on' ? [...multi.params.levels].sort((a, b) => a.radius - b.radius) : [{ radius: p.radius, bias: 1 }];
            const packed = new Float32Array(Math.max(1, levels.length) * 4);
            for (let i = 0; i < levels.length; i++) packed.set([Math.pow(2, levels[i].radius) * camera.scale, levels[i].bias, 0, 0], i * 4);
            const key = JSON.stringify(Array.from(packed));
            if (key !== this.levelsKey) { this.levels.destroy(); this.levels = this.context.createBuffer(packed, GPUBufferUsage.STORAGE); this.levelsKey = key; }
            settings.set([Math.pow(2, p.radius) * camera.scale, p.bias, p.samples, p.transparentThreshold], 52);
            settings.set([p.blurKernelSize, p.blurDepthBias, 0, 0], 56);
            settings.set([...Color.toRgbNormalized(p.color), includeOpaqueAo ? 1 : 0], 60);
            settings.set([levels.length, multi.name === 'on' ? 1 : 0, multi.name === 'on' ? multi.params.nearThreshold : 0, multi.name === 'on' ? multi.params.farThreshold : 0], 64);
        }
        if (props.sharpening.name === 'on') {
            const p = props.sharpening.params;
            settings.set([Math.pow(2, -(2 - 2 * Math.pow(p.sharpness, 0.25))), p.denoise ? 1 : 0, 0, 0], 68);
        }
        if (props.dof.name === 'on') {
            const p = props.dof.params;
            const center = p.center === 'scene-center' ? sceneCenter : camera.state.target;
            const viewCenter = Vec3.transformMat4(Vec3(), Vec3.scale(Vec3(), center, camera.scale), camera.view);
            settings.set([Math.max(1, Math.round(p.blurSize * pixelRatio)), p.blurSpread, (Vec3.distance(camera.state.position, center) + p.inFocus) * camera.scale, p.PPM * camera.scale], 72);
            settings.set([...viewCenter, p.mode === 'sphere' ? 1 : 0], 76);
        }
        if (props.shadow.name === 'on' && rendererProps) {
            const p = props.shadow.params;
            const light = getLight(rendererProps.light);
            const direction = Mat4.isZero(camera.headRotation) ? light.direction : getTransformedLightDirection(light, Mat4.invert(Mat4(), camera.headRotation));
            const packed = new Float32Array(Math.max(1, light.count) * 8);
            for (let i = 0; i < light.count; i++) { packed.set(direction.slice(i * 3, i * 3 + 3), i * 8); packed.set(light.color.slice(i * 3, i * 3 + 3), i * 8 + 4); }
            const key = JSON.stringify(Array.from(packed));
            if (key !== this.lightsKey) { this.lights.destroy(); this.lights = this.context.createBuffer(packed, GPUBufferUsage.STORAGE); this.lightsKey = key; }
            settings.set([p.maxDistance * camera.scale, p.tolerance * camera.scale, Math.max(1, p.steps), light.count], 80);
            settings.set([...Color.toRgbNormalized(rendererProps.ambientColor).map(c => c * rendererProps.ambientIntensity), camera.state.fog ? camera.fogNear : 0], 84);
            settings.set([...(transparentBackground ? [0, 0, 0] : Color.toRgbNormalized(rendererProps.backgroundColor)), camera.state.fog ? camera.fogFar : 0], 88);
        }
        const bloom = bloomEnabled && props.bloom.name === 'on' ? this.bloom.render(encoder, color, emissive, depth, picking, camera, props.bloom.params) : color;
        const aoDepth = aoEnabled && !reuseAo ? this.depthPyramid.render(encoder, depth, transparentDepth, aoWidth, aoHeight) : transparentDepth;
        let processedTransparentColor = transparentColor;
        const run = (name: string, pipeline: GPURenderPipeline, source: GPUTexture, target: GPUTexture, ao: GPUTexture, direction?: [number, number], outlines = color) => {
            let uniform = this.uniforms.get(name);
            if (!uniform) { uniform = this.context.device.createBuffer({ size: 384, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST }); this.uniforms.set(name, uniform); }
            settings[58] = direction?.[0] ?? 0; settings[59] = direction?.[1] ?? 0;
            this.context.device.queue.writeBuffer(uniform, 0, settings);
            const bindings = this.context.device.createBindGroup({ layout: this.layout, entries: [
                { binding: 0, resource: { buffer: uniform } }, { binding: 1, resource: source.createView() },
                { binding: 2, resource: this.sampler }, { binding: 3, resource: depth.createView() }, { binding: 4, resource: picking.createView() },
                { binding: 5, resource: ao.createView() }, { binding: 6, resource: { buffer: this.samples } }, { binding: 7, resource: { buffer: this.levels } }, { binding: 8, resource: outlines.createView() }, { binding: 9, resource: bloom.createView() }, { binding: 10, resource: { buffer: this.lights } },
                { binding: 11, resource: emissive.createView() },
                { binding: 12, resource: transparentDepth.createView() }, { binding: 13, resource: processedTransparentColor.createView() },
                { binding: 14, resource: aoDepth.createView() },
            ] });
            const pass = encoder.beginRenderPass({ label: `molstar-postprocess-${name}`, colorAttachments: [{ view: target.createView(), loadOp: 'clear', storeOp: 'store' }] });
            pass.setPipeline(pipeline); pass.setBindGroup(0, bindings); pass.draw(3); pass.end();
        };
        if (aoEnabled) {
            if (!reuseAo) {
                run('ssao', this.pipelines.ssao, color, this.aoTargets[0], color);
                run('ssaoBlurX', this.pipelines.ssaoBlur, color, this.aoTargets[1], this.aoTargets[0], [1, 0]);
                run('ssaoBlurY', this.pipelines.ssaoBlur, color, this.aoTargets[0], this.aoTargets[1], [0, 1]);
            }
            if (includeTransparentAo) run('ssaoTransparentColor', this.pipelines.ssaoTransparentColor, color, this.transparentAoColor!, this.aoTargets[0]);
        }
        let source = color, marked = false;
        for (const [i, stage] of stages.entries()) {
            if ((applyMarking || applyBackground) && !marked && (stage === 'dof' || stage === 'fxaa' || stage === 'smaa' || stage === 'sharpen')) { source = finish(source); marked = true; }
            if (stage === 'smaa' && props.antialiasing.name === 'smaa') { source = this.smaa.render(encoder, source, camera, props.antialiasing.params); continue; }
            if (stage === 'smaa') continue;
            const target = this.targets[i % 2];
            if (stage === 'outline') run('outlineEdges', this.pipelines.outlineEdges, source, this.outlineTarget!, color);
            run(stage, this.pipelines[stage], source, target, aoEnabled ? this.aoTargets[0] : color, undefined, stage === 'outline' ? this.outlineTarget! : color);
            source = target;
            if (stage === 'ssaoCompose' && includeTransparentAo) processedTransparentColor = this.transparentAoColor!;
        }
        return !marked ? finish(source) : source;
    }

    dispose() { this.depthPyramid.dispose(); this.smaa.dispose(); this.bloom.dispose(); this.outlineTarget?.destroy(); this.transparentAoColor?.destroy(); for (const u of this.uniforms.values()) u.destroy(); this.uniforms.clear(); this.samples.destroy(); this.levels.destroy(); this.lights.destroy(); for (const t of [...this.targets, ...this.aoTargets]) t.destroy(); this.targets.length = 0; this.aoTargets.length = 0; }
}
