/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { TextureSpec, Values, UniformSpec, DefineSpec } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Texture } from '../../mol-gl/webgl/texture';
import { ValueCell } from '../../mol-util';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { createComputeRenderable, ComputeRenderable } from '../../mol-gl/renderable';
import { Mat4, Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { DrawPass } from './draw';
import { ICamera } from '../../mol-canvas3d/camera';
import quad_vert from '../../mol-gl/shader/quad.vert';
import outlines_frag from '../../mol-gl/shader/outlines.frag';
import ssao_frag from '../../mol-gl/shader/ssao.frag';
import ssao_blur_frag from '../../mol-gl/shader/ssao-blur.frag';
import postprocessing_frag from '../../mol-gl/shader/postprocessing.frag';
import fxaa_frag from '../../mol-gl/shader/fxaa.frag';
import { Framebuffer } from '../../mol-gl/webgl/framebuffer';
import { Color } from '../../mol-util/color';

const OutlinesSchema = {
    ...QuadSchema,
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),

    uMaxPossibleViewZDiff: UniformSpec('f'),
};
type OutlinesRenderable = ComputeRenderable<Values<typeof OutlinesSchema>>

function getOutlinesRenderable(ctx: WebGLContext, depthTexture: Texture): OutlinesRenderable {
    const values: Values<typeof OutlinesSchema> = {
        ...QuadValues,
        tDepth: ValueCell.create(depthTexture),
        uTexSize: ValueCell.create(Vec2.create(depthTexture.getWidth(), depthTexture.getHeight())),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),

        uMaxPossibleViewZDiff: ValueCell.create(0.5),
    };

    const schema = { ...OutlinesSchema };
    const shaderCode = ShaderCode('outlines', quad_vert, outlines_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

const SsaoSchema = {
    ...QuadSchema,
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),

    uSamples: UniformSpec('v3[]'),
    dNSamples: DefineSpec('number'),

    uProjection: UniformSpec('m4'),
    uInvProjection: UniformSpec('m4'),

    uTexSize: UniformSpec('v2'),

    uRadius: UniformSpec('f'),
    uBias: UniformSpec('f'),
};

type SsaoRenderable = ComputeRenderable<Values<typeof SsaoSchema>>

function getSsaoRenderable(ctx: WebGLContext, depthTexture: Texture, nSamples: number): SsaoRenderable {
    const values: Values<typeof SsaoSchema> = {
        ...QuadValues,
        tDepth: ValueCell.create(depthTexture),

        uSamples: ValueCell.create(getSamples(nSamples)),
        dNSamples: ValueCell.create(nSamples),

        uProjection: ValueCell.create(Mat4.identity()),
        uInvProjection: ValueCell.create(Mat4.identity()),

        uTexSize: ValueCell.create(Vec2.create(ctx.gl.drawingBufferWidth, ctx.gl.drawingBufferHeight)),

        uRadius: ValueCell.create(8.0),
        uBias: ValueCell.create(0.025),
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

    uMaxPossibleViewZDiff: UniformSpec('f'),

    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    dOrthographic: DefineSpec('number'),
};

type SsaoBlurRenderable = ComputeRenderable<Values<typeof SsaoBlurSchema>>

function getSsaoBlurRenderable(ctx: WebGLContext, ssaoDepthTexture: Texture, blurKernelSize: number, direction: 'horizontal' | 'vertical'): SsaoBlurRenderable {
    const values: Values<typeof SsaoBlurSchema> = {
        ...QuadValues,
        tSsaoDepth: ValueCell.create(ssaoDepthTexture),
        uTexSize: ValueCell.create(Vec2.create(ssaoDepthTexture.getWidth(), ssaoDepthTexture.getHeight())),

        uKernel: ValueCell.create(getBlurKernel(blurKernelSize)),
        dOcclusionKernelSize: ValueCell.create(blurKernelSize),

        uBlurDirectionX: ValueCell.create(direction === 'horizontal' ? 1 : 0),
        uBlurDirectionY: ValueCell.create(direction === 'vertical' ? 1 : 0),

        uMaxPossibleViewZDiff: ValueCell.create(0.5),

        uNear: ValueCell.create(0.0),
        uFar: ValueCell.create(10000.0),
        dOrthographic: ValueCell.create(0),
    };

    const schema = { ...SsaoBlurSchema };
    const shaderCode = ShaderCode('ssao_blur', quad_vert, ssao_blur_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

function getBlurKernel(kernelSize: number): number[] {
    let sigma = kernelSize / 3.0;
    let halfKernelSize = Math.floor((kernelSize + 1) / 2);

    let kernel = [];
    for (let x = 0; x < halfKernelSize; x++) {
        kernel.push((1.0 / ((Math.sqrt(2 * Math.PI)) * sigma)) * Math.exp(-x * x / (2 * sigma * sigma)));
    }

    return kernel;
}

function getSamples(nSamples: number): number[] {
    let vectorSamples = [];
    for (let i = 0; i < nSamples; i++) {
        let v = Vec3();

        v[0] = Math.random() * 2.0 - 1.0;
        v[1] = Math.random() * 2.0 - 1.0;
        v[2] = Math.random();

        Vec3.normalize(v, v);

        Vec3.scale(v, v, Math.random());

        let scale = (i * i) / (nSamples * nSamples);
        scale = 0.1 + scale * (1.0 - 0.1);

        Vec3.scale(v, v, scale);

        vectorSamples.push(v);
    }

    let samples = [];
    for (let i = 0; i < nSamples; i++) {
        let v = vectorSamples[i];
        samples.push(v[0]);
        samples.push(v[1]);
        samples.push(v[2]);
    }

    return samples;
}

const PostprocessingSchema = {
    ...QuadSchema,
    tSsaoDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tOutlines: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    uFogNear: UniformSpec('f'),
    uFogFar: UniformSpec('f'),
    uFogColor: UniformSpec('v3'),

    uMaxPossibleViewZDiff: UniformSpec('f'),

    dOcclusionEnable: DefineSpec('boolean'),

    dOutlineEnable: DefineSpec('boolean'),
    uOutlineScale: UniformSpec('f'),
    uOutlineThreshold: UniformSpec('f'),

    dPackedDepth: DefineSpec('boolean'),
};
type PostprocessingRenderable = ComputeRenderable<Values<typeof PostprocessingSchema>>

function getPostprocessingRenderable(ctx: WebGLContext, colorTexture: Texture, depthTexture: Texture, packedDepth: boolean, outlinesTexture: Texture, ssaoDepthTexture: Texture): PostprocessingRenderable {
    const values: Values<typeof PostprocessingSchema> = {
        ...QuadValues,
        tSsaoDepth: ValueCell.create(ssaoDepthTexture),
        tColor: ValueCell.create(colorTexture),
        tDepth: ValueCell.create(depthTexture),
        tOutlines: ValueCell.create(outlinesTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),
        uFogNear: ValueCell.create(10000),
        uFogFar: ValueCell.create(10000),
        uFogColor: ValueCell.create(Vec3.create(1, 1, 1)),

        uMaxPossibleViewZDiff: ValueCell.create(0.5),

        dOcclusionEnable: ValueCell.create(false),

        dOutlineEnable: ValueCell.create(false),
        uOutlineScale: ValueCell.create(ctx.pixelRatio),
        uOutlineThreshold: ValueCell.create(0.8),

        dPackedDepth: ValueCell.create(packedDepth),
    };

    const schema = { ...PostprocessingSchema };
    const shaderCode = ShaderCode('postprocessing', quad_vert, postprocessing_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

export const PostprocessingParams = {
    occlusion: PD.MappedStatic('off', {
        on: PD.Group({
            samples: PD.Numeric(64, {min: 1, max: 256, step: 1}),
            radius: PD.Numeric(8.0, { min: 0.1, max: 64, step: 0.1 }),
            bias: PD.Numeric(0.025, { min: 0, max: 1, step: 0.001 }),
            kernelSize: PD.Numeric(13, { min: 1, max: 25, step: 2 }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Darken occluded crevices with the ambient occlusion effect' }),
    outline: PD.MappedStatic('off', {
        on: PD.Group({
            scale: PD.Numeric(1, { min: 1, max: 5, step: 1 }),
            threshold: PD.Numeric(0.1, { min: 0.01, max: 1, step: 0.01 }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Draw outline around 3D objects' }),
    antialiasing: PD.MappedStatic('on', {
        on: PD.Group({
            edgeThresholdMin:PD.Numeric(0.0312, { min: 0.0312, max: 0.0833, step: 0.0001 }, { description: 'Trims the algorithm from processing darks.' }),
            edgeThresholdMax: PD.Numeric(0.063, { min: 0.063, max: 0.333, step: 0.001 }, { description: 'The minimum amount of local contrast required to apply algorithm.' }),
            iterations: PD.Numeric(12, { min: 0, max: 32, step: 1 }, { description: 'Number of edge exploration steps.' }),
            subpixelQuality: PD.Numeric(1.00, { min: 0.00, max: 1.00, step: 0.01 }, { description: 'Choose the amount of sub-pixel aliasing removal.' }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Fast Approximate Anti-Aliasing (FXAA)' }),
};
export type PostprocessingProps = PD.Values<typeof PostprocessingParams>

export class PostprocessingPass {
    static isEnabled(props: PostprocessingProps) {
        return props.occlusion.name === 'on' || props.outline.name === 'on' || props.antialiasing.name === 'on';
    }

    readonly target: RenderTarget

    private readonly outlinesTarget: RenderTarget
    private readonly outlinesRenderable: OutlinesRenderable

    private readonly ssaoFramebuffer: Framebuffer
    private readonly ssaoBlurFirstPassFramebuffer: Framebuffer
    private readonly ssaoBlurSecondPassFramebuffer: Framebuffer

    private readonly ssaoDepthTexture: Texture
    private readonly ssaoDepthBlurProxyTexture: Texture

    private readonly ssaoRenderable: SsaoRenderable
    private readonly ssaoBlurFirstPassRenderable: SsaoBlurRenderable
    private readonly ssaoBlurSecondPassRenderable: SsaoBlurRenderable

    private nSamples: number
    private blurKernelSize: number

    private readonly renderable: PostprocessingRenderable

    constructor(private webgl: WebGLContext, drawPass: DrawPass) {
        const { colorTarget, depthTexture, packedDepth } = drawPass;
        const width = colorTarget.getWidth();
        const height = colorTarget.getHeight();

        this.nSamples = 64;
        this.blurKernelSize = 13;

        this.target = webgl.createRenderTarget(width, height, false, 'uint8', 'linear');

        this.outlinesTarget = webgl.createRenderTarget(width, height, false);
        this.outlinesRenderable = getOutlinesRenderable(webgl, depthTexture);

        this.ssaoFramebuffer = webgl.resources.framebuffer();
        this.ssaoBlurFirstPassFramebuffer = webgl.resources.framebuffer();
        this.ssaoBlurSecondPassFramebuffer = webgl.resources.framebuffer();

        this.ssaoDepthTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
        this.ssaoDepthTexture.define(width, height);
        this.ssaoDepthTexture.attachFramebuffer(this.ssaoFramebuffer, 'color0');

        this.ssaoDepthBlurProxyTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
        this.ssaoDepthBlurProxyTexture.define(width, height);
        this.ssaoDepthBlurProxyTexture.attachFramebuffer(this.ssaoBlurFirstPassFramebuffer, 'color0');

        this.ssaoDepthTexture.attachFramebuffer(this.ssaoBlurSecondPassFramebuffer, 'color0');

        this.ssaoRenderable = getSsaoRenderable(webgl, depthTexture, this.nSamples);
        this.ssaoBlurFirstPassRenderable = getSsaoBlurRenderable(webgl, this.ssaoDepthTexture, this.blurKernelSize, 'horizontal');
        this.ssaoBlurSecondPassRenderable = getSsaoBlurRenderable(webgl, this.ssaoDepthBlurProxyTexture, this.blurKernelSize, 'vertical');
        this.renderable = getPostprocessingRenderable(webgl, colorTarget.texture,  depthTexture, packedDepth, this.outlinesTarget.texture, this.ssaoDepthTexture);
    }

    setSize(width: number, height: number) {
        const [w, h] = this.renderable.values.uTexSize.ref.value;
        if (width !== w || height !== h) {
            this.target.setSize(width, height);
            this.outlinesTarget.setSize(width, height);
            this.ssaoDepthTexture.define(width, height);
            this.ssaoDepthBlurProxyTexture.define(width, height);

            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.outlinesRenderable.values.uTexSize, Vec2.set(this.outlinesRenderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.ssaoRenderable.values.uTexSize, Vec2.set(this.ssaoRenderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uTexSize, Vec2.set(this.ssaoRenderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uTexSize, Vec2.set(this.ssaoRenderable.values.uTexSize.ref.value, width, height));
        }
    }

    private updateState(camera: ICamera, backgroundColor: Color, props: PostprocessingProps) {
        let needsUpdateMain = false;
        let needsUpdateSsao = false;
        let needsUpdateSsaoBlur = false;

        let orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
        let outlinesEnabled = props.outline.name === 'on';
        let occlusionEnabled = props.occlusion.name === 'on';

        let invProjection = Mat4.identity();
        Mat4.invert(invProjection, camera.projection);

        if (props.occlusion.name === 'on') {
            ValueCell.updateIfChanged(this.ssaoRenderable.values.uProjection, camera.projection);
            ValueCell.updateIfChanged(this.ssaoRenderable.values.uInvProjection, invProjection);

            ValueCell.updateIfChanged(this.ssaoBlurFirstPassRenderable.values.uNear, camera.near);
            ValueCell.updateIfChanged(this.ssaoBlurSecondPassRenderable.values.uNear, camera.near);

            ValueCell.updateIfChanged(this.ssaoBlurFirstPassRenderable.values.uFar, camera.far);
            ValueCell.updateIfChanged(this.ssaoBlurSecondPassRenderable.values.uFar, camera.far);

            if (this.ssaoBlurFirstPassRenderable.values.dOrthographic.ref.value !== orthographic) { needsUpdateSsaoBlur = true; }
            ValueCell.updateIfChanged(this.ssaoBlurFirstPassRenderable.values.dOrthographic, orthographic);
            ValueCell.updateIfChanged(this.ssaoBlurSecondPassRenderable.values.dOrthographic, orthographic);

            if (this.nSamples !== props.occlusion.params.samples) {
                needsUpdateSsao = true;

                this.nSamples = props.occlusion.params.samples;
                ValueCell.updateIfChanged(this.ssaoRenderable.values.uSamples, getSamples(this.nSamples));
                ValueCell.updateIfChanged(this.ssaoRenderable.values.dNSamples, this.nSamples);
            }
            ValueCell.updateIfChanged(this.ssaoRenderable.values.uRadius, props.occlusion.params.radius);
            ValueCell.updateIfChanged(this.ssaoRenderable.values.uBias, props.occlusion.params.bias);

            if (this.blurKernelSize !== props.occlusion.params.kernelSize) {
                needsUpdateSsaoBlur = true;

                this.blurKernelSize = props.occlusion.params.kernelSize;
                let kernel = getBlurKernel(this.blurKernelSize);

                ValueCell.updateIfChanged(this.ssaoBlurFirstPassRenderable.values.uKernel, kernel);
                ValueCell.updateIfChanged(this.ssaoBlurSecondPassRenderable.values.uKernel, kernel);
                ValueCell.updateIfChanged(this.ssaoBlurFirstPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
                ValueCell.updateIfChanged(this.ssaoBlurSecondPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
            }

        }

        if (props.outline.name === 'on') {
            let maxPossibleViewZDiff = props.outline.params.threshold * (camera.fogFar - camera.near);

            ValueCell.updateIfChanged(this.outlinesRenderable.values.uNear, camera.near);
            ValueCell.updateIfChanged(this.outlinesRenderable.values.uFar, camera.far);
            ValueCell.updateIfChanged(this.outlinesRenderable.values.uMaxPossibleViewZDiff, maxPossibleViewZDiff);

            ValueCell.updateIfChanged(this.renderable.values.uMaxPossibleViewZDiff, maxPossibleViewZDiff);
            let fogColor = Vec3();
            Color.toVec3Normalized(fogColor, backgroundColor);
            ValueCell.updateIfChanged(this.renderable.values.uFogColor, fogColor);
            ValueCell.updateIfChanged(this.renderable.values.uOutlineScale, props.outline.params.scale - 1);
            ValueCell.updateIfChanged(this.renderable.values.uOutlineThreshold, props.outline.params.threshold);
        }

        ValueCell.updateIfChanged(this.renderable.values.uFar, camera.far);
        ValueCell.updateIfChanged(this.renderable.values.uNear, camera.near);
        ValueCell.updateIfChanged(this.renderable.values.uFogFar, camera.fogFar);
        ValueCell.updateIfChanged(this.renderable.values.uFogNear, camera.fogNear);
        if (this.renderable.values.dOrthographic.ref.value !== orthographic) { needsUpdateMain = true; }
        ValueCell.updateIfChanged(this.renderable.values.dOrthographic, orthographic);
        if (this.renderable.values.dOutlineEnable.ref.value !== outlinesEnabled) { needsUpdateMain = true; }
        ValueCell.updateIfChanged(this.renderable.values.dOutlineEnable, outlinesEnabled);
        if (this.renderable.values.dOcclusionEnable.ref.value !== occlusionEnabled) { needsUpdateMain = true; }
        ValueCell.updateIfChanged(this.renderable.values.dOcclusionEnable, occlusionEnabled);

        if (needsUpdateSsao) {
            this.ssaoRenderable.update();
        }

        if (needsUpdateSsaoBlur) {
            this.ssaoBlurFirstPassRenderable.update();
            this.ssaoBlurSecondPassRenderable.update();
        }

        if (needsUpdateMain) {
            this.renderable.update();
        }

        const { gl, state } = this.webgl;

        state.enable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);

        const { x, y, width, height } = camera.viewport;
        gl.viewport(x, y, width, height);
        gl.scissor(x, y, width, height);
    }

    render(camera: ICamera, toDrawingBuffer: boolean, backgroundColor: Color, props: PostprocessingProps) {
        this.updateState(camera, backgroundColor, props);

        if (props.outline.name === 'on') {
            this.outlinesTarget.bind();
            this.outlinesRenderable.render();
        }

        if (props.occlusion.name === 'on') {
            this.ssaoFramebuffer.bind();
            this.ssaoRenderable.render();

            this.ssaoBlurFirstPassFramebuffer.bind();
            this.ssaoBlurFirstPassRenderable.render();

            this.ssaoBlurSecondPassFramebuffer.bind();
            this.ssaoBlurSecondPassRenderable.render();
        }

        if (toDrawingBuffer) {
            this.webgl.unbindFramebuffer();
        } else {
            this.target.bind();
        }

        const { gl, state } = this.webgl;
        state.clearColor(0, 0, 0, 1);
        gl.clear(gl.COLOR_BUFFER_BIT);

        this.renderable.render();
    }
}

export class FxaaPass {
    static isEnabled(props: PostprocessingProps) {
        return props.antialiasing.name === 'on';
    }

    readonly target: RenderTarget
    private readonly renderable: FxaaRenderable

    constructor(private webgl: WebGLContext, private drawPass: DrawPass) {
        const { colorTarget } = drawPass;
        const width = colorTarget.getWidth();
        const height = colorTarget.getHeight();

        this.target = webgl.createRenderTarget(width, height, false);
        this.renderable = getFxaaRenderable(webgl, drawPass.colorTarget.texture);
    }

    setSize(width: number, height: number) {
        const [w, h] = [this.target.texture.getWidth(), this.target.texture.getHeight()];
        if (width !== w || height !== h) {
            this.target.setSize(width, height);
            ValueCell.update(this.renderable.values.uTexSizeInv, Vec2.set(this.renderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
        }
    }

    private updateState(camera: ICamera) {
        const { gl, state } = this.webgl;

        state.enable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);

        const { x, y, width, height } = camera.viewport;
        gl.viewport(x, y, width, height);
        gl.scissor(x, y, width, height);
    }

    render(camera: ICamera, toDrawingBuffer: boolean, props: PostprocessingProps) {
        if (props.antialiasing.name === 'off') return;

        const { values } = this.renderable;

        let needsUpdate = false;

        if (PostprocessingPass.isEnabled(props)) {
            ValueCell.update(this.renderable.values.tColor, this.drawPass.postprocessing.target.texture);
            needsUpdate = true;
        } else {
            ValueCell.update(this.renderable.values.tColor, this.drawPass.colorTarget.texture);
            needsUpdate = true;
        }

        const { edgeThresholdMin, edgeThresholdMax, iterations, subpixelQuality } = props.antialiasing.params;
        if (values.dEdgeThresholdMin.ref.value !== edgeThresholdMin) needsUpdate = true;
        ValueCell.updateIfChanged(values.dEdgeThresholdMin, edgeThresholdMin);
        if (values.dEdgeThresholdMax.ref.value !== edgeThresholdMax) needsUpdate = true;
        ValueCell.updateIfChanged(values.dEdgeThresholdMax, edgeThresholdMax);
        if (values.dIterations.ref.value !== iterations) needsUpdate = true;
        ValueCell.updateIfChanged(values.dIterations, iterations);
        if (values.dSubpixelQuality.ref.value !== subpixelQuality) needsUpdate = true;
        ValueCell.updateIfChanged(values.dSubpixelQuality, subpixelQuality);

        if (needsUpdate) {
            this.renderable.update();
        }

        if (toDrawingBuffer) {
            this.webgl.unbindFramebuffer();
        } else {
            this.target.bind();
        }

        this.updateState(camera);
        this.renderable.render();
    }
}

//

const FxaaSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSizeInv: UniformSpec('v2'),

    dEdgeThresholdMin: DefineSpec('number'),
    dEdgeThresholdMax: DefineSpec('number'),
    dIterations: DefineSpec('number'),
    dSubpixelQuality: DefineSpec('number'),
};
const FxaaShaderCode = ShaderCode('fxaa', quad_vert, fxaa_frag);
type FxaaRenderable = ComputeRenderable<Values<typeof FxaaSchema>>

function getFxaaRenderable(ctx: WebGLContext, colorTexture: Texture): FxaaRenderable {
    const width = colorTexture.getWidth();
    const height = colorTexture.getHeight();

    const values: Values<typeof FxaaSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),

        dEdgeThresholdMin: ValueCell.create(0.0312),
        dEdgeThresholdMax: ValueCell.create(0.125),
        dIterations: ValueCell.create(12),
        dSubpixelQuality: ValueCell.create(0.75),
    };

    const schema = { ...FxaaSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', FxaaShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}
