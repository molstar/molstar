/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { CopyRenderable, createCopyRenderable, QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { TextureSpec, Values, UniformSpec, DefineSpec } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Texture } from '../../mol-gl/webgl/texture';
import { deepEqual, ValueCell } from '../../mol-util';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { createComputeRenderable, ComputeRenderable } from '../../mol-gl/renderable';
import { Mat4, Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { DrawPass } from './draw';
import { ICamera } from '../../mol-canvas3d/camera';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { outlines_frag } from '../../mol-gl/shader/outlines.frag';
import { ssao_frag } from '../../mol-gl/shader/ssao.frag';
import { ssaoBlur_frag } from '../../mol-gl/shader/ssao-blur.frag';
import { postprocessing_frag } from '../../mol-gl/shader/postprocessing.frag';
import { Framebuffer } from '../../mol-gl/webgl/framebuffer';
import { Color } from '../../mol-util/color';
import { FxaaParams, FxaaPass } from './fxaa';
import { SmaaParams, SmaaPass } from './smaa';
import { isTimingMode } from '../../mol-util/debug';
import { BackgroundParams, BackgroundPass } from './background';
import { AssetManager } from '../../mol-util/assets';
import { Light } from '../../mol-gl/renderer';
import { shadows_frag } from '../../mol-gl/shader/shadows.frag';
import { CasParams, CasPass } from './cas';

const OutlinesSchema = {
    ...QuadSchema,
    tDepthOpaque: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    uInvProjection: UniformSpec('m4'),

    uOutlineThreshold: UniformSpec('f'),
    dTransparentOutline: DefineSpec('boolean'),
};
type OutlinesRenderable = ComputeRenderable<Values<typeof OutlinesSchema>>

function getOutlinesRenderable(ctx: WebGLContext, depthTextureOpaque: Texture, depthTextureTransparent: Texture, transparentOutline: boolean): OutlinesRenderable {
    const width = depthTextureOpaque.getWidth();
    const height = depthTextureOpaque.getHeight();

    const values: Values<typeof OutlinesSchema> = {
        ...QuadValues,
        tDepthOpaque: ValueCell.create(depthTextureOpaque),
        tDepthTransparent: ValueCell.create(depthTextureTransparent),
        uTexSize: ValueCell.create(Vec2.create(width, height)),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),
        uInvProjection: ValueCell.create(Mat4.identity()),

        uOutlineThreshold: ValueCell.create(0.33),
        dTransparentOutline: ValueCell.create(transparentOutline),
    };

    const schema = { ...OutlinesSchema };
    const shaderCode = ShaderCode('outlines', quad_vert, outlines_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

const ShadowsSchema = {
    ...QuadSchema,
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    uProjection: UniformSpec('m4'),
    uInvProjection: UniformSpec('m4'),
    uBounds: UniformSpec('v4'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),

    dSteps: DefineSpec('number'),
    uMaxDistance: UniformSpec('f'),
    uTolerance: UniformSpec('f'),
    uBias: UniformSpec('f'),

    uLightDirection: UniformSpec('v3[]'),
    uLightColor: UniformSpec('v3[]'),
    dLightCount: DefineSpec('number'),
};
type ShadowsRenderable = ComputeRenderable<Values<typeof ShadowsSchema>>

function getShadowsRenderable(ctx: WebGLContext, depthTexture: Texture): ShadowsRenderable {
    const width = depthTexture.getWidth();
    const height = depthTexture.getHeight();

    const values: Values<typeof ShadowsSchema> = {
        ...QuadValues,
        tDepth: ValueCell.create(depthTexture),
        uTexSize: ValueCell.create(Vec2.create(width, height)),

        uProjection: ValueCell.create(Mat4.identity()),
        uInvProjection: ValueCell.create(Mat4.identity()),
        uBounds: ValueCell.create(Vec4()),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),

        dSteps: ValueCell.create(1),
        uMaxDistance: ValueCell.create(3.0),
        uTolerance: ValueCell.create(1.0),
        uBias: ValueCell.create(0.6),

        uLightDirection: ValueCell.create([]),
        uLightColor: ValueCell.create([]),
        dLightCount: ValueCell.create(0),
    };

    const schema = { ...ShadowsSchema };
    const shaderCode = ShaderCode('shadows', quad_vert, shadows_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

const SsaoSchema = {
    ...QuadSchema,
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepthHalf: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepthQuarter: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),

    uSamples: UniformSpec('v3[]'),
    dNSamples: DefineSpec('number'),

    uProjection: UniformSpec('m4'),
    uInvProjection: UniformSpec('m4'),
    uBounds: UniformSpec('v4'),

    uTexSize: UniformSpec('v2'),

    uRadius: UniformSpec('f'),
    uBias: UniformSpec('f'),

    dMultiScale: DefineSpec('boolean'),
    dLevels: DefineSpec('number'),
    uLevelRadius: UniformSpec('f[]'),
    uLevelBias: UniformSpec('f[]'),
    uNearThreshold: UniformSpec('f'),
    uFarThreshold: UniformSpec('f'),
};

type SsaoRenderable = ComputeRenderable<Values<typeof SsaoSchema>>

function getSsaoRenderable(ctx: WebGLContext, depthTexture: Texture, depthHalfTexture: Texture, depthQuarterTexture: Texture): SsaoRenderable {
    const values: Values<typeof SsaoSchema> = {
        ...QuadValues,
        tDepth: ValueCell.create(depthTexture),
        tDepthHalf: ValueCell.create(depthHalfTexture),
        tDepthQuarter: ValueCell.create(depthQuarterTexture),

        uSamples: ValueCell.create(getSamples(32)),
        dNSamples: ValueCell.create(32),

        uProjection: ValueCell.create(Mat4.identity()),
        uInvProjection: ValueCell.create(Mat4.identity()),
        uBounds: ValueCell.create(Vec4()),

        uTexSize: ValueCell.create(Vec2.create(ctx.gl.drawingBufferWidth, ctx.gl.drawingBufferHeight)),

        uRadius: ValueCell.create(Math.pow(2, 5)),
        uBias: ValueCell.create(0.8),

        dMultiScale: ValueCell.create(false),
        dLevels: ValueCell.create(3),
        uLevelRadius: ValueCell.create([Math.pow(2, 2), Math.pow(2, 5), Math.pow(2, 8)]),
        uLevelBias: ValueCell.create([0.8, 0.8, 0.8]),
        uNearThreshold: ValueCell.create(10.0),
        uFarThreshold: ValueCell.create(1500.0),
    };

    const schema = { ...SsaoSchema };
    const shaderCode = ShaderCode('ssao', quad_vert, ssao_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

const SsaoBlurSchema = {
    ...QuadSchema,
    tSsaoDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    uKernel: UniformSpec('f[]'),
    dOcclusionKernelSize: DefineSpec('number'),

    uBlurDirectionX: UniformSpec('f'),
    uBlurDirectionY: UniformSpec('f'),

    uInvProjection: UniformSpec('m4'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    uBounds: UniformSpec('v4'),
    dOrthographic: DefineSpec('number'),
};

type SsaoBlurRenderable = ComputeRenderable<Values<typeof SsaoBlurSchema>>

function getSsaoBlurRenderable(ctx: WebGLContext, ssaoDepthTexture: Texture, direction: 'horizontal' | 'vertical'): SsaoBlurRenderable {
    const values: Values<typeof SsaoBlurSchema> = {
        ...QuadValues,
        tSsaoDepth: ValueCell.create(ssaoDepthTexture),
        uTexSize: ValueCell.create(Vec2.create(ssaoDepthTexture.getWidth(), ssaoDepthTexture.getHeight())),

        uKernel: ValueCell.create(getBlurKernel(15)),
        dOcclusionKernelSize: ValueCell.create(15),

        uBlurDirectionX: ValueCell.create(direction === 'horizontal' ? 1 : 0),
        uBlurDirectionY: ValueCell.create(direction === 'vertical' ? 1 : 0),

        uInvProjection: ValueCell.create(Mat4.identity()),
        uNear: ValueCell.create(0.0),
        uFar: ValueCell.create(10000.0),
        uBounds: ValueCell.create(Vec4()),
        dOrthographic: ValueCell.create(0),
    };

    const schema = { ...SsaoBlurSchema };
    const shaderCode = ShaderCode('ssao_blur', quad_vert, ssaoBlur_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

function getBlurKernel(kernelSize: number): number[] {
    const sigma = kernelSize / 3.0;
    const halfKernelSize = Math.floor((kernelSize + 1) / 2);

    const kernel = [];
    for (let x = 0; x < halfKernelSize; x++) {
        kernel.push((1.0 / ((Math.sqrt(2 * Math.PI)) * sigma)) * Math.exp(-x * x / (2 * sigma * sigma)));
    }

    return kernel;
}

const RandomHemisphereVector: Vec3[] = [];
for (let i = 0; i < 256; i++) {
    const v = Vec3();
    v[0] = Math.random() * 2.0 - 1.0;
    v[1] = Math.random() * 2.0 - 1.0;
    v[2] = Math.random();
    Vec3.normalize(v, v);
    Vec3.scale(v, v, Math.random());
    RandomHemisphereVector.push(v);
}

function getSamples(nSamples: number): number[] {
    const samples = [];
    for (let i = 0; i < nSamples; i++) {
        let scale = (i * i + 2.0 * i + 1) / (nSamples * nSamples);
        scale = 0.1 + scale * (1.0 - 0.1);

        samples.push(RandomHemisphereVector[i][0] * scale);
        samples.push(RandomHemisphereVector[i][1] * scale);
        samples.push(RandomHemisphereVector[i][2] * scale);
    }

    return samples;
}

const PostprocessingSchema = {
    ...QuadSchema,
    tSsaoDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepthOpaque: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tShadows: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tOutlines: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    uFogNear: UniformSpec('f'),
    uFogFar: UniformSpec('f'),
    uFogColor: UniformSpec('v3'),
    uOutlineColor: UniformSpec('v3'),
    uOcclusionColor: UniformSpec('v3'),
    uTransparentBackground: UniformSpec('b'),

    dOcclusionEnable: DefineSpec('boolean'),
    uOcclusionOffset: UniformSpec('v2'),

    dShadowEnable: DefineSpec('boolean'),

    dOutlineEnable: DefineSpec('boolean'),
    dOutlineScale: DefineSpec('number'),
    dTransparentOutline: DefineSpec('boolean'),
};
type PostprocessingRenderable = ComputeRenderable<Values<typeof PostprocessingSchema>>

function getPostprocessingRenderable(ctx: WebGLContext, colorTexture: Texture, depthTextureOpaque: Texture, depthTextureTransparent: Texture, shadowsTexture: Texture, outlinesTexture: Texture, ssaoDepthTexture: Texture, transparentOutline: boolean): PostprocessingRenderable {
    const values: Values<typeof PostprocessingSchema> = {
        ...QuadValues,
        tSsaoDepth: ValueCell.create(ssaoDepthTexture),
        tColor: ValueCell.create(colorTexture),
        tDepthOpaque: ValueCell.create(depthTextureOpaque),
        tDepthTransparent: ValueCell.create(depthTextureTransparent),
        tShadows: ValueCell.create(shadowsTexture),
        tOutlines: ValueCell.create(outlinesTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),
        uFogNear: ValueCell.create(10000),
        uFogFar: ValueCell.create(10000),
        uFogColor: ValueCell.create(Vec3.create(1, 1, 1)),
        uOutlineColor: ValueCell.create(Vec3.create(0, 0, 0)),
        uOcclusionColor: ValueCell.create(Vec3.create(0, 0, 0)),
        uTransparentBackground: ValueCell.create(false),

        dOcclusionEnable: ValueCell.create(true),
        uOcclusionOffset: ValueCell.create(Vec2.create(0, 0)),

        dShadowEnable: ValueCell.create(false),

        dOutlineEnable: ValueCell.create(false),
        dOutlineScale: ValueCell.create(1),
        dTransparentOutline: ValueCell.create(transparentOutline),
    };

    const schema = { ...PostprocessingSchema };
    const shaderCode = ShaderCode('postprocessing', quad_vert, postprocessing_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

export const PostprocessingParams = {
    occlusion: PD.MappedStatic('on', {
        on: PD.Group({
            samples: PD.Numeric(32, { min: 1, max: 256, step: 1 }),
            multiScale: PD.MappedStatic('off', {
                on: PD.Group({
                    levels: PD.ObjectList({
                        radius: PD.Numeric(5, { min: 0, max: 20, step: 0.1 }, { description: 'Final occlusion radius is 2^x' }),
                        bias: PD.Numeric(1, { min: 0, max: 3, step: 0.1 }),
                    }, o => `${o.radius}, ${o.bias}`, { defaultValue: [
                        { radius: 2, bias: 1 },
                        { radius: 5, bias: 1 },
                        { radius: 8, bias: 1 },
                        { radius: 11, bias: 1 },
                    ] }),
                    nearThreshold: PD.Numeric(10, { min: 0, max: 50, step: 1 }),
                    farThreshold: PD.Numeric(1500, { min: 0, max: 10000, step: 100 }),
                }),
                off: PD.Group({})
            }, { cycle: true }),
            radius: PD.Numeric(5, { min: 0, max: 20, step: 0.1 }, { description: 'Final occlusion radius is 2^x', hideIf: p => p?.multiScale.name === 'on' }),
            bias: PD.Numeric(0.8, { min: 0, max: 3, step: 0.1 }),
            blurKernelSize: PD.Numeric(15, { min: 1, max: 25, step: 2 }),
            resolutionScale: PD.Numeric(1, { min: 0.1, max: 1, step: 0.05 }, { description: 'Adjust resolution of occlusion calculation' }),
            color: PD.Color(Color(0x000000)),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Darken occluded crevices with the ambient occlusion effect' }),
    shadow: PD.MappedStatic('off', {
        on: PD.Group({
            steps: PD.Numeric(1, { min: 1, max: 64, step: 1 }),
            bias: PD.Numeric(0.6, { min: 0.0, max: 1.0, step: 0.01 }),
            maxDistance: PD.Numeric(3, { min: 0, max: 256, step: 1 }),
            tolerance: PD.Numeric(1.0, { min: 0.0, max: 10.0, step: 0.1 }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Simplistic shadows' }),
    outline: PD.MappedStatic('off', {
        on: PD.Group({
            scale: PD.Numeric(1, { min: 1, max: 5, step: 1 }),
            threshold: PD.Numeric(0.33, { min: 0.01, max: 1, step: 0.01 }),
            color: PD.Color(Color(0x000000)),
            includeTransparent: PD.Boolean(true, { description: 'Whether to show outline for transparent objects' }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Draw outline around 3D objects' }),
    antialiasing: PD.MappedStatic('smaa', {
        fxaa: PD.Group(FxaaParams),
        smaa: PD.Group(SmaaParams),
        off: PD.Group({})
    }, { options: [['fxaa', 'FXAA'], ['smaa', 'SMAA'], ['off', 'Off']], description: 'Smooth pixel edges' }),
    sharpening: PD.MappedStatic('off', {
        on: PD.Group(CasParams),
        off: PD.Group({})
    }, { cycle: true, description: 'Contrast Adaptive Sharpening' }),
    background: PD.Group(BackgroundParams, { isFlat: true }),
};

export type PostprocessingProps = PD.Values<typeof PostprocessingParams>

type Levels = {
    count: number
    radius: number[]
    bias: number[]
}

function getLevels(props: { radius: number, bias: number }[], levels?: Levels): Levels {
    const count = props.length;
    const { radius, bias } = levels || {
        radius: (new Array(count * 3)).fill(0),
        bias: (new Array(count * 3)).fill(0),
    };
    props = props.slice().sort((a, b) => a.radius - b.radius);
    for (let i = 0; i < count; ++i) {
        const p = props[i];
        radius[i] = Math.pow(2, p.radius);
        bias[i] = p.bias;
    }
    return { count, radius, bias };
}

export class PostprocessingPass {
    static isEnabled(props: PostprocessingProps) {
        return props.occlusion.name === 'on' || props.shadow.name === 'on' || props.outline.name === 'on' || props.background.variant.name !== 'off';
    }

    static isTransparentOutlineEnabled(props: PostprocessingProps) {
        return props.outline.name === 'on' && props.outline.params.includeTransparent;
    }

    readonly target: RenderTarget;

    private readonly outlinesTarget: RenderTarget;
    private readonly outlinesRenderable: OutlinesRenderable;

    private readonly shadowsTarget: RenderTarget;
    private readonly shadowsRenderable: ShadowsRenderable;

    private readonly ssaoFramebuffer: Framebuffer;
    private readonly ssaoBlurFirstPassFramebuffer: Framebuffer;
    private readonly ssaoBlurSecondPassFramebuffer: Framebuffer;

    private readonly downsampledDepthTarget: RenderTarget;
    private readonly downsampleDepthRenderable: CopyRenderable;

    private readonly depthHalfTarget: RenderTarget;
    private readonly depthHalfRenderable: CopyRenderable;

    private readonly depthQuarterTarget: RenderTarget;
    private readonly depthQuarterRenderable: CopyRenderable;

    private readonly ssaoDepthTexture: Texture;
    private readonly ssaoDepthBlurProxyTexture: Texture;

    private readonly ssaoRenderable: SsaoRenderable;
    private readonly ssaoBlurFirstPassRenderable: SsaoBlurRenderable;
    private readonly ssaoBlurSecondPassRenderable: SsaoBlurRenderable;

    private nSamples: number;
    private blurKernelSize: number;

    private readonly renderable: PostprocessingRenderable;

    private ssaoScale: number;
    private calcSsaoScale(resolutionScale: number) {
        // downscale ssao for high pixel-ratios
        return Math.min(1, 1 / this.webgl.pixelRatio) * resolutionScale;
    }

    private levels: { radius: number, bias: number }[];

    private readonly bgColor = Vec3();
    readonly background: BackgroundPass;

    constructor(private readonly webgl: WebGLContext, assetManager: AssetManager, private readonly drawPass: DrawPass) {
        const { colorTarget, depthTextureTransparent, depthTextureOpaque } = drawPass;
        const width = colorTarget.getWidth();
        const height = colorTarget.getHeight();

        this.nSamples = 1;
        this.blurKernelSize = 1;
        this.ssaoScale = this.calcSsaoScale(1);
        this.levels = [];

        // needs to be linear for anti-aliasing pass
        this.target = webgl.createRenderTarget(width, height, false, 'uint8', 'linear');

        this.outlinesTarget = webgl.createRenderTarget(width, height, false);
        this.outlinesRenderable = getOutlinesRenderable(webgl, depthTextureOpaque, depthTextureTransparent, true);

        this.shadowsTarget = webgl.createRenderTarget(width, height, false);
        this.shadowsRenderable = getShadowsRenderable(webgl, depthTextureOpaque);

        this.ssaoFramebuffer = webgl.resources.framebuffer();
        this.ssaoBlurFirstPassFramebuffer = webgl.resources.framebuffer();
        this.ssaoBlurSecondPassFramebuffer = webgl.resources.framebuffer();

        const sw = Math.floor(width * this.ssaoScale);
        const sh = Math.floor(height * this.ssaoScale);

        const hw = Math.max(1, Math.floor(sw * 0.5));
        const hh = Math.max(1, Math.floor(sh * 0.5));

        const qw = Math.max(1, Math.floor(sw * 0.25));
        const qh = Math.max(1, Math.floor(sh * 0.25));

        this.downsampledDepthTarget = webgl.createRenderTarget(sw, sh, false, 'float32', 'linear', webgl.isWebGL2 ? 'alpha' : 'rgba');
        this.downsampleDepthRenderable = createCopyRenderable(webgl, depthTextureOpaque);

        const depthTexture = this.ssaoScale === 1 ? depthTextureOpaque : this.downsampledDepthTarget.texture;

        this.depthHalfTarget = webgl.createRenderTarget(hw, hh, false, 'float32', 'linear', webgl.isWebGL2 ? 'alpha' : 'rgba');
        this.depthHalfRenderable = createCopyRenderable(webgl, depthTexture);

        this.depthQuarterTarget = webgl.createRenderTarget(qw, qh, false, 'float32', 'linear', webgl.isWebGL2 ? 'alpha' : 'rgba');
        this.depthQuarterRenderable = createCopyRenderable(webgl, this.depthHalfTarget.texture);

        this.ssaoDepthTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
        this.ssaoDepthTexture.define(sw, sh);
        this.ssaoDepthTexture.attachFramebuffer(this.ssaoFramebuffer, 'color0');

        this.ssaoDepthBlurProxyTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
        this.ssaoDepthBlurProxyTexture.define(sw, sh);
        this.ssaoDepthBlurProxyTexture.attachFramebuffer(this.ssaoBlurFirstPassFramebuffer, 'color0');

        this.ssaoDepthTexture.attachFramebuffer(this.ssaoBlurSecondPassFramebuffer, 'color0');

        this.ssaoRenderable = getSsaoRenderable(webgl, depthTexture, this.depthHalfTarget.texture, this.depthQuarterTarget.texture);
        this.ssaoBlurFirstPassRenderable = getSsaoBlurRenderable(webgl, this.ssaoDepthTexture, 'horizontal');
        this.ssaoBlurSecondPassRenderable = getSsaoBlurRenderable(webgl, this.ssaoDepthBlurProxyTexture, 'vertical');
        this.renderable = getPostprocessingRenderable(webgl, colorTarget.texture, depthTextureOpaque, depthTextureTransparent, this.shadowsTarget.texture, this.outlinesTarget.texture, this.ssaoDepthTexture, true);

        this.background = new BackgroundPass(webgl, assetManager, width, height);
    }

    setSize(width: number, height: number) {
        const [w, h] = this.renderable.values.uTexSize.ref.value;
        const ssaoScale = this.calcSsaoScale(1);

        if (width !== w || height !== h || this.ssaoScale !== ssaoScale) {
            this.ssaoScale = ssaoScale;

            this.target.setSize(width, height);
            this.outlinesTarget.setSize(width, height);
            this.shadowsTarget.setSize(width, height);

            const sw = Math.floor(width * this.ssaoScale);
            const sh = Math.floor(height * this.ssaoScale);
            this.downsampledDepthTarget.setSize(sw, sh);
            this.ssaoDepthTexture.define(sw, sh);
            this.ssaoDepthBlurProxyTexture.define(sw, sh);

            const hw = Math.max(1, Math.floor(sw * 0.5));
            const hh = Math.max(1, Math.floor(sh * 0.5));
            this.depthHalfTarget.setSize(hw, hh);

            const qw = Math.max(1, Math.floor(sw * 0.25));
            const qh = Math.max(1, Math.floor(sh * 0.25));
            this.depthQuarterTarget.setSize(qw, qh);

            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.outlinesRenderable.values.uTexSize, Vec2.set(this.outlinesRenderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.shadowsRenderable.values.uTexSize, Vec2.set(this.shadowsRenderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.downsampleDepthRenderable.values.uTexSize, Vec2.set(this.downsampleDepthRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.depthHalfRenderable.values.uTexSize, Vec2.set(this.depthHalfRenderable.values.uTexSize.ref.value, hw, hh));
            ValueCell.update(this.depthQuarterRenderable.values.uTexSize, Vec2.set(this.depthQuarterRenderable.values.uTexSize.ref.value, qw, qh));
            ValueCell.update(this.ssaoRenderable.values.uTexSize, Vec2.set(this.ssaoRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uTexSize, Vec2.set(this.ssaoBlurFirstPassRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uTexSize, Vec2.set(this.ssaoBlurSecondPassRenderable.values.uTexSize.ref.value, sw, sh));

            const depthTexture = this.ssaoScale === 1 ? this.drawPass.depthTextureOpaque : this.downsampledDepthTarget.texture;
            ValueCell.update(this.depthHalfRenderable.values.tColor, depthTexture);
            ValueCell.update(this.ssaoRenderable.values.tDepth, depthTexture);

            this.depthHalfRenderable.update();
            this.ssaoRenderable.update();

            this.background.setSize(width, height);
        }
    }

    private updateState(camera: ICamera, transparentBackground: boolean, backgroundColor: Color, props: PostprocessingProps, light: Light) {
        let needsUpdateShadows = false;
        let needsUpdateMain = false;
        let needsUpdateSsao = false;
        let needsUpdateSsaoBlur = false;
        let needsUpdateDepthHalf = false;
        let needsUpdateOutlines = false;

        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
        const outlinesEnabled = props.outline.name === 'on';
        const shadowsEnabled = props.shadow.name === 'on';
        const occlusionEnabled = props.occlusion.name === 'on';

        const invProjection = Mat4.identity();
        Mat4.invert(invProjection, camera.projection);

        const [w, h] = this.renderable.values.uTexSize.ref.value;
        const v = camera.viewport;

        if (props.occlusion.name === 'on') {
            ValueCell.update(this.ssaoRenderable.values.uProjection, camera.projection);
            ValueCell.update(this.ssaoRenderable.values.uInvProjection, invProjection);

            const b = this.ssaoRenderable.values.uBounds;
            const s = this.ssaoScale;
            Vec4.set(b.ref.value,
                Math.floor(v.x * s) / (w * s),
                Math.floor(v.y * s) / (h * s),
                Math.ceil((v.x + v.width) * s) / (w * s),
                Math.ceil((v.y + v.height) * s) / (h * s)
            );
            ValueCell.update(b, b.ref.value);
            ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uBounds, b.ref.value);
            ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uBounds, b.ref.value);

            ValueCell.updateIfChanged(this.ssaoBlurFirstPassRenderable.values.uNear, camera.near);
            ValueCell.updateIfChanged(this.ssaoBlurSecondPassRenderable.values.uNear, camera.near);

            ValueCell.updateIfChanged(this.ssaoBlurFirstPassRenderable.values.uFar, camera.far);
            ValueCell.updateIfChanged(this.ssaoBlurSecondPassRenderable.values.uFar, camera.far);

            ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uInvProjection, invProjection);
            ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uInvProjection, invProjection);

            if (this.ssaoBlurFirstPassRenderable.values.dOrthographic.ref.value !== orthographic) {
                needsUpdateSsaoBlur = true;
                ValueCell.update(this.ssaoBlurFirstPassRenderable.values.dOrthographic, orthographic);
                ValueCell.update(this.ssaoBlurSecondPassRenderable.values.dOrthographic, orthographic);
            }

            if (this.nSamples !== props.occlusion.params.samples) {
                needsUpdateSsao = true;

                this.nSamples = props.occlusion.params.samples;
                ValueCell.update(this.ssaoRenderable.values.uSamples, getSamples(this.nSamples));
                ValueCell.updateIfChanged(this.ssaoRenderable.values.dNSamples, this.nSamples);
            }

            const multiScale = props.occlusion.params.multiScale.name === 'on';
            if (this.ssaoRenderable.values.dMultiScale.ref.value !== multiScale) {
                needsUpdateSsao = true;
                ValueCell.update(this.ssaoRenderable.values.dMultiScale, multiScale);
            }

            if (props.occlusion.params.multiScale.name === 'on') {
                const mp = props.occlusion.params.multiScale.params;
                if (!deepEqual(this.levels, mp.levels)) {
                    needsUpdateSsao = true;

                    this.levels = mp.levels;
                    const levels = getLevels(mp.levels);
                    ValueCell.updateIfChanged(this.ssaoRenderable.values.dLevels, levels.count);

                    ValueCell.update(this.ssaoRenderable.values.uLevelRadius, levels.radius);
                    ValueCell.update(this.ssaoRenderable.values.uLevelBias, levels.bias);
                }
                ValueCell.updateIfChanged(this.ssaoRenderable.values.uNearThreshold, mp.nearThreshold);
                ValueCell.updateIfChanged(this.ssaoRenderable.values.uFarThreshold, mp.farThreshold);
            } else {
                ValueCell.updateIfChanged(this.ssaoRenderable.values.uRadius, Math.pow(2, props.occlusion.params.radius));
            }
            ValueCell.updateIfChanged(this.ssaoRenderable.values.uBias, props.occlusion.params.bias);

            if (this.blurKernelSize !== props.occlusion.params.blurKernelSize) {
                needsUpdateSsaoBlur = true;

                this.blurKernelSize = props.occlusion.params.blurKernelSize;
                const kernel = getBlurKernel(this.blurKernelSize);

                ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uKernel, kernel);
                ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uKernel, kernel);
                ValueCell.update(this.ssaoBlurFirstPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
                ValueCell.update(this.ssaoBlurSecondPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
            }

            const ssaoScale = this.calcSsaoScale(props.occlusion.params.resolutionScale);
            if (this.ssaoScale !== ssaoScale) {
                needsUpdateSsao = true;
                needsUpdateDepthHalf = true;

                this.ssaoScale = ssaoScale;

                const sw = Math.floor(w * this.ssaoScale);
                const sh = Math.floor(h * this.ssaoScale);
                this.downsampledDepthTarget.setSize(sw, sh);
                this.ssaoDepthTexture.define(sw, sh);
                this.ssaoDepthBlurProxyTexture.define(sw, sh);

                const hw = Math.floor(sw * 0.5);
                const hh = Math.floor(sh * 0.5);
                this.depthHalfTarget.setSize(hw, hh);

                const qw = Math.floor(sw * 0.25);
                const qh = Math.floor(sh * 0.25);
                this.depthQuarterTarget.setSize(qw, qh);

                const depthTexture = this.ssaoScale === 1 ? this.drawPass.depthTextureOpaque : this.downsampledDepthTarget.texture;
                ValueCell.update(this.depthHalfRenderable.values.tColor, depthTexture);
                ValueCell.update(this.ssaoRenderable.values.tDepth, depthTexture);

                ValueCell.update(this.ssaoRenderable.values.tDepthHalf, this.depthHalfTarget.texture);
                ValueCell.update(this.ssaoRenderable.values.tDepthQuarter, this.depthQuarterTarget.texture);

                ValueCell.update(this.downsampleDepthRenderable.values.uTexSize, Vec2.set(this.downsampleDepthRenderable.values.uTexSize.ref.value, sw, sh));
                ValueCell.update(this.depthHalfRenderable.values.uTexSize, Vec2.set(this.depthHalfRenderable.values.uTexSize.ref.value, hw, hh));
                ValueCell.update(this.depthQuarterRenderable.values.uTexSize, Vec2.set(this.depthQuarterRenderable.values.uTexSize.ref.value, qw, qh));
                ValueCell.update(this.ssaoRenderable.values.uTexSize, Vec2.set(this.ssaoRenderable.values.uTexSize.ref.value, sw, sh));
                ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uTexSize, Vec2.set(this.ssaoBlurFirstPassRenderable.values.uTexSize.ref.value, sw, sh));
                ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uTexSize, Vec2.set(this.ssaoBlurSecondPassRenderable.values.uTexSize.ref.value, sw, sh));
            }

            ValueCell.update(this.renderable.values.uOcclusionColor, Color.toVec3Normalized(this.renderable.values.uOcclusionColor.ref.value, props.occlusion.params.color));
        }

        if (props.shadow.name === 'on') {
            ValueCell.update(this.shadowsRenderable.values.uProjection, camera.projection);
            ValueCell.update(this.shadowsRenderable.values.uInvProjection, invProjection);

            Vec4.set(this.shadowsRenderable.values.uBounds.ref.value,
                v.x / w,
                v.y / h,
                (v.x + v.width) / w,
                (v.y + v.height) / h
            );
            ValueCell.update(this.shadowsRenderable.values.uBounds, this.shadowsRenderable.values.uBounds.ref.value);

            ValueCell.updateIfChanged(this.shadowsRenderable.values.uNear, camera.near);
            ValueCell.updateIfChanged(this.shadowsRenderable.values.uFar, camera.far);
            if (this.shadowsRenderable.values.dOrthographic.ref.value !== orthographic) {
                ValueCell.update(this.shadowsRenderable.values.dOrthographic, orthographic);
                needsUpdateShadows = true;
            }

            ValueCell.updateIfChanged(this.shadowsRenderable.values.uMaxDistance, props.shadow.params.maxDistance);
            ValueCell.updateIfChanged(this.shadowsRenderable.values.uTolerance, props.shadow.params.tolerance);
            ValueCell.updateIfChanged(this.shadowsRenderable.values.uBias, props.shadow.params.bias);
            if (this.shadowsRenderable.values.dSteps.ref.value !== props.shadow.params.steps) {
                ValueCell.update(this.shadowsRenderable.values.dSteps, props.shadow.params.steps);
                needsUpdateShadows = true;
            }

            ValueCell.update(this.shadowsRenderable.values.uLightDirection, light.direction);
            ValueCell.update(this.shadowsRenderable.values.uLightColor, light.color);
            if (this.shadowsRenderable.values.dLightCount.ref.value !== light.count) {
                ValueCell.update(this.shadowsRenderable.values.dLightCount, light.count);
                needsUpdateShadows = true;
            }
        }

        if (props.outline.name === 'on') {
            const transparentOutline = props.outline.params.includeTransparent ?? true;
            const outlineScale = Math.max(1, Math.round(props.outline.params.scale * this.webgl.pixelRatio)) - 1;
            const outlineThreshold = 50 * props.outline.params.threshold * this.webgl.pixelRatio;

            ValueCell.updateIfChanged(this.outlinesRenderable.values.uNear, camera.near);
            ValueCell.updateIfChanged(this.outlinesRenderable.values.uFar, camera.far);
            ValueCell.update(this.outlinesRenderable.values.uInvProjection, invProjection);
            if (this.outlinesRenderable.values.dTransparentOutline.ref.value !== transparentOutline) {
                needsUpdateOutlines = true;
                ValueCell.update(this.outlinesRenderable.values.dTransparentOutline, transparentOutline);
            }
            if (this.outlinesRenderable.values.dOrthographic.ref.value !== orthographic) {
                needsUpdateOutlines = true;
                ValueCell.update(this.outlinesRenderable.values.dOrthographic, orthographic);
            }
            ValueCell.updateIfChanged(this.outlinesRenderable.values.uOutlineThreshold, outlineThreshold);

            ValueCell.update(this.renderable.values.uOutlineColor, Color.toVec3Normalized(this.renderable.values.uOutlineColor.ref.value, props.outline.params.color));

            if (this.renderable.values.dOutlineScale.ref.value !== outlineScale) {
                needsUpdateMain = true;
                ValueCell.update(this.renderable.values.dOutlineScale, outlineScale);
            }
            if (this.renderable.values.dTransparentOutline.ref.value !== transparentOutline) {
                needsUpdateMain = true;
                ValueCell.update(this.renderable.values.dTransparentOutline, transparentOutline);
            }
        }

        ValueCell.updateIfChanged(this.renderable.values.uFar, camera.far);
        ValueCell.updateIfChanged(this.renderable.values.uNear, camera.near);
        ValueCell.updateIfChanged(this.renderable.values.uFogFar, camera.fogFar);
        ValueCell.updateIfChanged(this.renderable.values.uFogNear, camera.fogNear);
        ValueCell.update(this.renderable.values.uFogColor, Color.toVec3Normalized(this.renderable.values.uFogColor.ref.value, backgroundColor));
        ValueCell.updateIfChanged(this.renderable.values.uTransparentBackground, transparentBackground);
        if (this.renderable.values.dOrthographic.ref.value !== orthographic) {
            needsUpdateMain = true;
            ValueCell.update(this.renderable.values.dOrthographic, orthographic);
        }

        if (this.renderable.values.dOutlineEnable.ref.value !== outlinesEnabled) {
            needsUpdateMain = true;
            ValueCell.update(this.renderable.values.dOutlineEnable, outlinesEnabled);
        }
        if (this.renderable.values.dShadowEnable.ref.value !== shadowsEnabled) {
            needsUpdateMain = true;
            ValueCell.update(this.renderable.values.dShadowEnable, shadowsEnabled);
        }
        if (this.renderable.values.dOcclusionEnable.ref.value !== occlusionEnabled) {
            needsUpdateMain = true;
            ValueCell.update(this.renderable.values.dOcclusionEnable, occlusionEnabled);
        }

        if (needsUpdateOutlines) {
            this.outlinesRenderable.update();
        }

        if (needsUpdateShadows) {
            this.shadowsRenderable.update();
        }

        if (needsUpdateSsao) {
            this.ssaoRenderable.update();
        }

        if (needsUpdateSsaoBlur) {
            this.ssaoBlurFirstPassRenderable.update();
            this.ssaoBlurSecondPassRenderable.update();
        }

        if (needsUpdateDepthHalf) {
            this.depthHalfRenderable.update();
        }

        if (needsUpdateMain) {
            this.renderable.update();
        }

        const { gl, state } = this.webgl;

        state.enable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);
    }

    private occlusionOffset: [x: number, y: number] = [0, 0];
    setOcclusionOffset(x: number, y: number) {
        this.occlusionOffset[0] = x;
        this.occlusionOffset[1] = y;
        ValueCell.update(this.renderable.values.uOcclusionOffset, Vec2.set(this.renderable.values.uOcclusionOffset.ref.value, x, y));
    }

    private transparentBackground = false;
    setTransparentBackground(value: boolean) {
        this.transparentBackground = value;
    }

    render(camera: ICamera, toDrawingBuffer: boolean, transparentBackground: boolean, backgroundColor: Color, props: PostprocessingProps, light: Light) {
        if (isTimingMode) this.webgl.timer.mark('PostprocessingPass.render');
        this.updateState(camera, transparentBackground, backgroundColor, props, light);

        const { gl, state } = this.webgl;
        const { x, y, width, height } = camera.viewport;

        // don't render occlusion if offset is given,
        // which will reuse the existing occlusion
        if (props.occlusion.name === 'on' && this.occlusionOffset[0] === 0 && this.occlusionOffset[1] === 0) {
            if (isTimingMode) this.webgl.timer.mark('SSAO.render');
            const sx = Math.floor(x * this.ssaoScale);
            const sy = Math.floor(y * this.ssaoScale);
            const sw = Math.ceil(width * this.ssaoScale);
            const sh = Math.ceil(height * this.ssaoScale);

            state.viewport(sx, sy, sw, sh);
            state.scissor(sx, sy, sw, sh);

            if (this.ssaoScale < 1) {
                if (isTimingMode) this.webgl.timer.mark('SSAO.downsample');
                this.downsampledDepthTarget.bind();
                this.downsampleDepthRenderable.render();
                if (isTimingMode) this.webgl.timer.markEnd('SSAO.downsample');
            }

            if (isTimingMode) this.webgl.timer.mark('SSAO.half');
            this.depthHalfTarget.bind();
            this.depthHalfRenderable.render();
            if (isTimingMode) this.webgl.timer.markEnd('SSAO.half');

            if (isTimingMode) this.webgl.timer.mark('SSAO.quarter');
            this.depthQuarterTarget.bind();
            this.depthQuarterRenderable.render();
            if (isTimingMode) this.webgl.timer.markEnd('SSAO.quarter');

            this.ssaoFramebuffer.bind();
            this.ssaoRenderable.render();

            this.ssaoBlurFirstPassFramebuffer.bind();
            this.ssaoBlurFirstPassRenderable.render();

            this.ssaoBlurSecondPassFramebuffer.bind();
            this.ssaoBlurSecondPassRenderable.render();
            if (isTimingMode) this.webgl.timer.markEnd('SSAO.render');
        }

        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        if (props.outline.name === 'on') {
            this.outlinesTarget.bind();
            this.outlinesRenderable.render();
        }

        if (props.shadow.name === 'on') {
            this.shadowsTarget.bind();
            this.shadowsRenderable.render();
        }

        if (toDrawingBuffer) {
            this.webgl.unbindFramebuffer();
        } else {
            this.target.bind();
        }

        this.background.update(camera, props.background);
        if (this.background.isEnabled(props.background)) {
            if (this.transparentBackground) {
                state.clearColor(0, 0, 0, 0);
            } else {
                Color.toVec3Normalized(this.bgColor, backgroundColor);
                state.clearColor(this.bgColor[0], this.bgColor[1], this.bgColor[2], 1);
            }
            gl.clear(gl.COLOR_BUFFER_BIT);
            state.enable(gl.BLEND);
            state.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
            this.background.render();
        } else {
            state.clearColor(0, 0, 0, 1);
            gl.clear(gl.COLOR_BUFFER_BIT);
        }

        this.renderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('PostprocessingPass.render');
    }
}

export class AntialiasingPass {
    static isEnabled(props: PostprocessingProps) {
        return props.antialiasing.name !== 'off';
    }

    readonly target: RenderTarget;
    private readonly internalTarget: RenderTarget;

    private readonly fxaa: FxaaPass;
    private readonly smaa: SmaaPass;
    private readonly cas: CasPass;

    constructor(webgl: WebGLContext, private drawPass: DrawPass) {
        const { colorTarget } = drawPass;
        const width = colorTarget.getWidth();
        const height = colorTarget.getHeight();

        this.target = webgl.createRenderTarget(width, height, false);
        this.internalTarget = webgl.createRenderTarget(width, height, false);

        this.fxaa = new FxaaPass(webgl, this.target.texture);
        this.smaa = new SmaaPass(webgl, this.target.texture);
        this.cas = new CasPass(webgl, this.target.texture);
    }

    setSize(width: number, height: number) {
        const w = this.target.texture.getWidth();
        const h = this.target.texture.getHeight();

        if (width !== w || height !== h) {
            this.target.setSize(width, height);
            this.internalTarget.setSize(width, height);
            this.fxaa.setSize(width, height);
            if (this.smaa.supported) this.smaa.setSize(width, height);
            this.cas.setSize(width, height);
        }
    }

    private _renderFxaa(camera: ICamera, target: RenderTarget | undefined, props: PostprocessingProps) {
        if (props.antialiasing.name !== 'fxaa') return;

        const input = PostprocessingPass.isEnabled(props)
            ? this.drawPass.postprocessing.target.texture
            : this.drawPass.colorTarget.texture;
        this.fxaa.update(input, props.antialiasing.params);
        this.fxaa.render(camera.viewport, target);
    }

    private _renderSmaa(camera: ICamera, target: RenderTarget | undefined, props: PostprocessingProps) {
        if (props.antialiasing.name !== 'smaa') return;

        const input = PostprocessingPass.isEnabled(props)
            ? this.drawPass.postprocessing.target.texture
            : this.drawPass.colorTarget.texture;
        this.smaa.update(input, props.antialiasing.params);
        this.smaa.render(camera.viewport, target);
    }

    private _renderAntialiasing(camera: ICamera, target: RenderTarget | undefined, props: PostprocessingProps) {
        if (props.antialiasing.name === 'fxaa') {
            this._renderFxaa(camera, target, props);
        } else if (props.antialiasing.name === 'smaa') {
            this._renderSmaa(camera, target, props);
        }
    }

    private _renderCas(camera: ICamera, target: RenderTarget | undefined, props: PostprocessingProps) {
        if (props.sharpening.name !== 'on') return;

        const input = props.antialiasing.name !== 'off'
            ? this.internalTarget.texture
            : PostprocessingPass.isEnabled(props)
                ? this.drawPass.postprocessing.target.texture
                : this.drawPass.colorTarget.texture;
        this.cas.update(input, props.sharpening.params);
        this.cas.render(camera.viewport, target);
    }

    render(camera: ICamera, toDrawingBuffer: boolean, props: PostprocessingProps) {
        if (props.antialiasing.name === 'off' && props.sharpening.name === 'off') return;

        if (props.antialiasing.name === 'smaa' && !this.smaa.supported) {
            console.error('SMAA not supported, missing "HTMLImageElement"');
            return;
        }

        const target = toDrawingBuffer ? undefined : this.target;
        if (props.sharpening.name === 'off') {
            this._renderAntialiasing(camera, target, props);
        } else if (props.antialiasing.name === 'off') {
            this._renderCas(camera, target, props);
        } else {
            this._renderAntialiasing(camera, this.internalTarget, props);
            this._renderCas(camera, target, props);
        }
    }
}
