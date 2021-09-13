/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { TextureSpec, Values, UniformSpec, DefineSpec } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Texture } from '../../mol-gl/webgl/texture';
import { ValueCell } from '../../mol-util';
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
    uBounds: UniformSpec('v4'),

    uTexSize: UniformSpec('v2'),

    uRadius: UniformSpec('f'),
    uBias: UniformSpec('f'),
};

type SsaoRenderable = ComputeRenderable<Values<typeof SsaoSchema>>

function getSsaoRenderable(ctx: WebGLContext, depthTexture: Texture): SsaoRenderable {
    const values: Values<typeof SsaoSchema> = {
        ...QuadValues,
        tDepth: ValueCell.create(depthTexture),

        uSamples: ValueCell.create([0.0, 0.0, 1.0]),
        dNSamples: ValueCell.create(1),

        uProjection: ValueCell.create(Mat4.identity()),
        uInvProjection: ValueCell.create(Mat4.identity()),
        uBounds: ValueCell.create(Vec4()),

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
    uBounds: UniformSpec('v4'),
    dOrthographic: DefineSpec('number'),
};

type SsaoBlurRenderable = ComputeRenderable<Values<typeof SsaoBlurSchema>>

function getSsaoBlurRenderable(ctx: WebGLContext, ssaoDepthTexture: Texture, direction: 'horizontal' | 'vertical'): SsaoBlurRenderable {
    const values: Values<typeof SsaoBlurSchema> = {
        ...QuadValues,
        tSsaoDepth: ValueCell.create(ssaoDepthTexture),
        uTexSize: ValueCell.create(Vec2.create(ssaoDepthTexture.getWidth(), ssaoDepthTexture.getHeight())),

        uKernel: ValueCell.create([0.0]),
        dOcclusionKernelSize: ValueCell.create(1),

        uBlurDirectionX: ValueCell.create(direction === 'horizontal' ? 1 : 0),
        uBlurDirectionY: ValueCell.create(direction === 'vertical' ? 1 : 0),

        uMaxPossibleViewZDiff: ValueCell.create(0.5),

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

function getSamples(vectorSamples: Vec3[], nSamples: number): number[] {
    const samples = [];
    for (let i = 0; i < nSamples; i++) {
        let scale = (i * i + 2.0 * i + 1) / (nSamples * nSamples);
        scale = 0.1 + scale * (1.0 - 0.1);

        samples.push(vectorSamples[i][0] * scale);
        samples.push(vectorSamples[i][1] * scale);
        samples.push(vectorSamples[i][2] * scale);
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
    uTransparentBackground: UniformSpec('b'),

    uMaxPossibleViewZDiff: UniformSpec('f'),

    dOcclusionEnable: DefineSpec('boolean'),

    dOutlineEnable: DefineSpec('boolean'),
    dOutlineScale: DefineSpec('number'),
    uOutlineThreshold: UniformSpec('f'),
};
type PostprocessingRenderable = ComputeRenderable<Values<typeof PostprocessingSchema>>

function getPostprocessingRenderable(ctx: WebGLContext, colorTexture: Texture, depthTexture: Texture, outlinesTexture: Texture, ssaoDepthTexture: Texture): PostprocessingRenderable {
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
        uTransparentBackground: ValueCell.create(false),

        uMaxPossibleViewZDiff: ValueCell.create(0.5),

        dOcclusionEnable: ValueCell.create(false),

        dOutlineEnable: ValueCell.create(false),
        dOutlineScale: ValueCell.create(1),
        uOutlineThreshold: ValueCell.create(0.33),
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
            radius: PD.Numeric(5, { min: 0, max: 10, step: 0.1 }, { description: 'Final radius is 2^x.' }),
            bias: PD.Numeric(0.8, { min: 0, max: 3, step: 0.1 }),
            blurKernelSize: PD.Numeric(15, { min: 1, max: 25, step: 2 }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Darken occluded crevices with the ambient occlusion effect' }),
    outline: PD.MappedStatic('off', {
        on: PD.Group({
            scale: PD.Numeric(1, { min: 1, max: 5, step: 1 }),
            threshold: PD.Numeric(0.33, { min: 0.01, max: 1, step: 0.01 }),
        }),
        off: PD.Group({})
    }, { cycle: true, description: 'Draw outline around 3D objects' }),
    antialiasing: PD.MappedStatic('smaa', {
        fxaa: PD.Group(FxaaParams),
        smaa: PD.Group(SmaaParams),
        off: PD.Group({})
    }, { options: [['fxaa', 'FXAA'], ['smaa', 'SMAA'], ['off', 'Off']], description: 'Smooth pixel edges' }),
};
export type PostprocessingProps = PD.Values<typeof PostprocessingParams>

export class PostprocessingPass {
    static isEnabled(props: PostprocessingProps) {
        return props.occlusion.name === 'on' || props.outline.name === 'on';
    }

    readonly target: RenderTarget

    private readonly outlinesTarget: RenderTarget
    private readonly outlinesRenderable: OutlinesRenderable

    private readonly randomHemisphereVector: Vec3[]
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

    private ssaoScale: number
    private calcSsaoScale() {
        // downscale ssao for high pixel-ratios
        return Math.min(1, 1 / this.webgl.pixelRatio);
    }

    constructor(private webgl: WebGLContext, drawPass: DrawPass) {
        this.ssaoScale = this.calcSsaoScale();

        const { colorTarget, depthTexture } = drawPass;
        const width = colorTarget.getWidth();
        const height = colorTarget.getHeight();

        this.nSamples = 1;
        this.blurKernelSize = 1;

        // needs to be linear for anti-aliasing pass
        this.target = webgl.createRenderTarget(width, height, false, 'uint8', 'linear');

        this.outlinesTarget = webgl.createRenderTarget(width, height, false);
        this.outlinesRenderable = getOutlinesRenderable(webgl, depthTexture);

        this.randomHemisphereVector = [];
        for (let i = 0; i < 256; i++) {
            const v = Vec3();
            v[0] = Math.random() * 2.0 - 1.0;
            v[1] = Math.random() * 2.0 - 1.0;
            v[2] = Math.random();
            Vec3.normalize(v, v);
            Vec3.scale(v, v, Math.random());
            this.randomHemisphereVector.push(v);
        }
        this.ssaoFramebuffer = webgl.resources.framebuffer();
        this.ssaoBlurFirstPassFramebuffer = webgl.resources.framebuffer();
        this.ssaoBlurSecondPassFramebuffer = webgl.resources.framebuffer();

        const sw = Math.floor(width * this.ssaoScale);
        const sh = Math.floor(height * this.ssaoScale);

        this.ssaoDepthTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
        this.ssaoDepthTexture.define(sw, sh);
        this.ssaoDepthTexture.attachFramebuffer(this.ssaoFramebuffer, 'color0');

        this.ssaoDepthBlurProxyTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
        this.ssaoDepthBlurProxyTexture.define(sw, sh);
        this.ssaoDepthBlurProxyTexture.attachFramebuffer(this.ssaoBlurFirstPassFramebuffer, 'color0');

        this.ssaoDepthTexture.attachFramebuffer(this.ssaoBlurSecondPassFramebuffer, 'color0');

        this.ssaoRenderable = getSsaoRenderable(webgl, depthTexture);
        this.ssaoBlurFirstPassRenderable = getSsaoBlurRenderable(webgl, this.ssaoDepthTexture, 'horizontal');
        this.ssaoBlurSecondPassRenderable = getSsaoBlurRenderable(webgl, this.ssaoDepthBlurProxyTexture, 'vertical');
        this.renderable = getPostprocessingRenderable(webgl, colorTarget.texture, depthTexture, this.outlinesTarget.texture, this.ssaoDepthTexture);
    }

    setSize(width: number, height: number) {
        const [w, h] = this.renderable.values.uTexSize.ref.value;
        const ssaoScale = this.calcSsaoScale();

        if (width !== w || height !== h || this.ssaoScale !== ssaoScale) {
            this.ssaoScale = ssaoScale;

            const sw = Math.floor(width * this.ssaoScale);
            const sh = Math.floor(height * this.ssaoScale);
            this.target.setSize(width, height);
            this.outlinesTarget.setSize(width, height);
            this.ssaoDepthTexture.define(sw, sh);
            this.ssaoDepthBlurProxyTexture.define(sw, sh);

            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.outlinesRenderable.values.uTexSize, Vec2.set(this.outlinesRenderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.ssaoRenderable.values.uTexSize, Vec2.set(this.ssaoRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uTexSize, Vec2.set(this.ssaoBlurFirstPassRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uTexSize, Vec2.set(this.ssaoBlurSecondPassRenderable.values.uTexSize.ref.value, sw, sh));
        }
    }

    private updateState(camera: ICamera, transparentBackground: boolean, backgroundColor: Color, props: PostprocessingProps) {
        let needsUpdateMain = false;
        let needsUpdateSsao = false;
        let needsUpdateSsaoBlur = false;

        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
        const outlinesEnabled = props.outline.name === 'on';
        const occlusionEnabled = props.occlusion.name === 'on';

        const invProjection = Mat4.identity();
        Mat4.invert(invProjection, camera.projection);

        if (props.occlusion.name === 'on') {
            ValueCell.update(this.ssaoRenderable.values.uProjection, camera.projection);
            ValueCell.update(this.ssaoRenderable.values.uInvProjection, invProjection);

            const [w, h] = this.renderable.values.uTexSize.ref.value;
            const b = this.ssaoRenderable.values.uBounds;
            const v = camera.viewport;
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

            if (this.ssaoBlurFirstPassRenderable.values.dOrthographic.ref.value !== orthographic) {
                needsUpdateSsaoBlur = true;
            }
            ValueCell.updateIfChanged(this.ssaoBlurFirstPassRenderable.values.dOrthographic, orthographic);
            ValueCell.updateIfChanged(this.ssaoBlurSecondPassRenderable.values.dOrthographic, orthographic);

            if (this.nSamples !== props.occlusion.params.samples) {
                needsUpdateSsao = true;

                this.nSamples = props.occlusion.params.samples;
                ValueCell.update(this.ssaoRenderable.values.uSamples, getSamples(this.randomHemisphereVector, this.nSamples));
                ValueCell.updateIfChanged(this.ssaoRenderable.values.dNSamples, this.nSamples);
            }
            ValueCell.updateIfChanged(this.ssaoRenderable.values.uRadius, Math.pow(2, props.occlusion.params.radius));
            ValueCell.updateIfChanged(this.ssaoRenderable.values.uBias, props.occlusion.params.bias);

            if (this.blurKernelSize !== props.occlusion.params.blurKernelSize) {
                needsUpdateSsaoBlur = true;

                this.blurKernelSize = props.occlusion.params.blurKernelSize;
                const kernel = getBlurKernel(this.blurKernelSize);

                ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uKernel, kernel);
                ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uKernel, kernel);
                ValueCell.updateIfChanged(this.ssaoBlurFirstPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
                ValueCell.updateIfChanged(this.ssaoBlurSecondPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
            }

        }

        if (props.outline.name === 'on') {
            let { threshold } = props.outline.params;
            // orthographic needs lower threshold
            if (camera.state.mode === 'orthographic') threshold /= 5;
            const factor = Math.pow(1000, threshold) / 1000;
            // use radiusMax for stable outlines when zooming
            const maxPossibleViewZDiff = factor * camera.state.radiusMax;
            const outlineScale = props.outline.params.scale - 1;

            ValueCell.updateIfChanged(this.outlinesRenderable.values.uNear, camera.near);
            ValueCell.updateIfChanged(this.outlinesRenderable.values.uFar, camera.far);
            ValueCell.updateIfChanged(this.outlinesRenderable.values.uMaxPossibleViewZDiff, maxPossibleViewZDiff);

            ValueCell.updateIfChanged(this.renderable.values.uMaxPossibleViewZDiff, maxPossibleViewZDiff);
            if (this.renderable.values.dOutlineScale.ref.value !== outlineScale) { needsUpdateMain = true; }
            ValueCell.updateIfChanged(this.renderable.values.dOutlineScale, outlineScale);
        }

        ValueCell.updateIfChanged(this.renderable.values.uFar, camera.far);
        ValueCell.updateIfChanged(this.renderable.values.uNear, camera.near);
        ValueCell.updateIfChanged(this.renderable.values.uFogFar, camera.fogFar);
        ValueCell.updateIfChanged(this.renderable.values.uFogNear, camera.fogNear);
        ValueCell.update(this.renderable.values.uFogColor, Color.toVec3Normalized(this.renderable.values.uFogColor.ref.value, backgroundColor));
        ValueCell.updateIfChanged(this.renderable.values.uTransparentBackground, transparentBackground);
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

    render(camera: ICamera, toDrawingBuffer: boolean, transparentBackground: boolean, backgroundColor: Color, props: PostprocessingProps) {
        this.updateState(camera, transparentBackground, backgroundColor, props);

        if (props.outline.name === 'on') {
            this.outlinesTarget.bind();
            this.outlinesRenderable.render();
        }

        if (props.occlusion.name === 'on') {
            const { x, y, width, height } = camera.viewport;
            const sx = Math.floor(x * this.ssaoScale);
            const sy = Math.floor(y * this.ssaoScale);
            const sw = Math.ceil(width * this.ssaoScale);
            const sh = Math.ceil(height * this.ssaoScale);
            this.webgl.gl.viewport(sx, sy, sw, sh);
            this.webgl.gl.scissor(sx, sy, sw, sh);

            this.ssaoFramebuffer.bind();
            this.ssaoRenderable.render();

            this.ssaoBlurFirstPassFramebuffer.bind();
            this.ssaoBlurFirstPassRenderable.render();

            this.ssaoBlurSecondPassFramebuffer.bind();
            this.ssaoBlurSecondPassRenderable.render();

            this.webgl.gl.viewport(x, y, width, height);
            this.webgl.gl.scissor(x, y, width, height);
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

export class AntialiasingPass {
    static isEnabled(props: PostprocessingProps) {
        return props.antialiasing.name !== 'off';
    }

    readonly target: RenderTarget
    private readonly fxaa: FxaaPass
    private readonly smaa: SmaaPass

    constructor(webgl: WebGLContext, private drawPass: DrawPass) {
        const { colorTarget } = drawPass;
        const width = colorTarget.getWidth();
        const height = colorTarget.getHeight();

        this.target = webgl.createRenderTarget(width, height, false);
        this.fxaa = new FxaaPass(webgl, this.target.texture);
        this.smaa = new SmaaPass(webgl, this.target.texture);
    }

    setSize(width: number, height: number) {
        const w = this.target.texture.getWidth();
        const h = this.target.texture.getHeight();

        if (width !== w || height !== h) {
            this.target.setSize(width, height);
            this.fxaa.setSize(width, height);
            if (this.smaa.supported) this.smaa.setSize(width, height);
        }
    }

    private _renderFxaa(camera: ICamera, toDrawingBuffer: boolean, props: PostprocessingProps) {
        if (props.antialiasing.name !== 'fxaa') return;

        const input = PostprocessingPass.isEnabled(props)
            ? this.drawPass.postprocessing.target.texture
            : this.drawPass.colorTarget.texture;
        this.fxaa.update(input, props.antialiasing.params);
        this.fxaa.render(camera.viewport, toDrawingBuffer ? undefined : this.target);
    }

    private _renderSmaa(camera: ICamera, toDrawingBuffer: boolean, props: PostprocessingProps) {
        if (props.antialiasing.name !== 'smaa') return;

        const input = PostprocessingPass.isEnabled(props)
            ? this.drawPass.postprocessing.target.texture
            : this.drawPass.colorTarget.texture;
        this.smaa.update(input, props.antialiasing.params);
        this.smaa.render(camera.viewport, toDrawingBuffer ? undefined : this.target);
    }

    render(camera: ICamera, toDrawingBuffer: boolean, props: PostprocessingProps) {
        if (props.antialiasing.name === 'off') return;

        if (props.antialiasing.name === 'fxaa') {
            this._renderFxaa(camera, toDrawingBuffer, props);
        } else if (props.antialiasing.name === 'smaa') {
            if (!this.smaa.supported) {
                throw new Error('SMAA not supported, missing "HTMLImageElement"');
            }
            this._renderSmaa(camera, toDrawingBuffer, props);
        }
    }
}

