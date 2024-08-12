/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
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
import { ICamera } from '../../mol-canvas3d/camera';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { ssao_frag } from '../../mol-gl/shader/ssao.frag';
import { ssaoBlur_frag } from '../../mol-gl/shader/ssao-blur.frag';
import { Framebuffer } from '../../mol-gl/webgl/framebuffer';
import { Color } from '../../mol-util/color';
import { isTimingMode } from '../../mol-util/debug';
import { PostprocessingProps } from './postprocessing';

export const SsaoParams = {
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
    blurDepthBias: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }),
    resolutionScale: PD.Numeric(1, { min: 0.1, max: 1, step: 0.05 }, { description: 'Adjust resolution of occlusion calculation' }),
    color: PD.Color(Color(0x000000)),
};

export type SsaoProps = PD.Values<typeof SsaoParams>

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

export class SsaoPass {
    static isEnabled(props: PostprocessingProps) {
        return props.occlusion.name !== 'off';
    }

    readonly target: RenderTarget;

    private readonly framebuffer: Framebuffer;
    private readonly blurFirstPassFramebuffer: Framebuffer;
    private readonly blurSecondPassFramebuffer: Framebuffer;

    private readonly downsampledDepthTarget: RenderTarget;
    private readonly downsampleDepthRenderable: CopyRenderable;

    private readonly depthHalfTarget: RenderTarget;
    private readonly depthHalfRenderable: CopyRenderable;

    private readonly depthQuarterTarget: RenderTarget;
    private readonly depthQuarterRenderable: CopyRenderable;

    readonly ssaoDepthTexture: Texture;
    private readonly depthBlurProxyTexture: Texture;

    private readonly renderable: SsaoRenderable;
    private readonly blurFirstPassRenderable: SsaoBlurRenderable;
    private readonly blurSecondPassRenderable: SsaoBlurRenderable;

    private depthTexture: Texture;

    private nSamples: number;
    private blurKernelSize: number;
    private texSize: [number, number];

    private ssaoScale: number;
    private calcSsaoScale(resolutionScale: number) {
        // downscale ssao for high pixel-ratios
        return Math.min(1, 1 / this.webgl.pixelRatio) * resolutionScale;
    }

    private levels: { radius: number, bias: number }[];

    private getDepthTexture() {
        return this.ssaoScale === 1 ? this.depthTexture : this.downsampledDepthTarget.texture;
    }

    constructor(private readonly webgl: WebGLContext, width: number, height: number, packedDepth: boolean, depthTexture: Texture) {
        const { textureFloatLinear } = webgl.extensions;

        this.depthTexture = depthTexture;

        this.nSamples = 1;
        this.blurKernelSize = 1;
        this.ssaoScale = this.calcSsaoScale(1);
        this.texSize = [width, height];
        this.levels = [];

        this.framebuffer = webgl.resources.framebuffer();
        this.blurFirstPassFramebuffer = webgl.resources.framebuffer();
        this.blurSecondPassFramebuffer = webgl.resources.framebuffer();

        const sw = Math.floor(width * this.ssaoScale);
        const sh = Math.floor(height * this.ssaoScale);

        const hw = Math.max(1, Math.floor(sw * 0.5));
        const hh = Math.max(1, Math.floor(sh * 0.5));

        const qw = Math.max(1, Math.floor(sw * 0.25));
        const qh = Math.max(1, Math.floor(sh * 0.25));

        const filter = textureFloatLinear ? 'linear' : 'nearest';

        this.downsampledDepthTarget = packedDepth
            ? webgl.createRenderTarget(sw, sh, false, 'uint8', filter, 'rgba')
            : webgl.createRenderTarget(sw, sh, false, 'float32', filter, webgl.isWebGL2 ? 'alpha' : 'rgba');
        this.downsampleDepthRenderable = createCopyRenderable(webgl, depthTexture);

        this.depthHalfTarget = packedDepth
            ? webgl.createRenderTarget(hw, hh, false, 'uint8', filter, 'rgba')
            : webgl.createRenderTarget(hw, hh, false, 'float32', filter, webgl.isWebGL2 ? 'alpha' : 'rgba');
        this.depthHalfRenderable = createCopyRenderable(webgl, this.getDepthTexture());

        this.depthQuarterTarget = packedDepth
            ? webgl.createRenderTarget(qw, qh, false, 'uint8', filter, 'rgba')
            : webgl.createRenderTarget(qw, qh, false, 'float32', filter, webgl.isWebGL2 ? 'alpha' : 'rgba');
        this.depthQuarterRenderable = createCopyRenderable(webgl, this.depthHalfTarget.texture);

        this.ssaoDepthTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', filter);
        this.ssaoDepthTexture.define(sw, sh);
        this.ssaoDepthTexture.attachFramebuffer(this.framebuffer, 'color0');

        this.depthBlurProxyTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', filter);
        this.depthBlurProxyTexture.define(sw, sh);
        this.depthBlurProxyTexture.attachFramebuffer(this.blurFirstPassFramebuffer, 'color0');

        this.ssaoDepthTexture.attachFramebuffer(this.blurSecondPassFramebuffer, 'color0');

        this.renderable = getSsaoRenderable(webgl, this.getDepthTexture(), this.depthHalfTarget.texture, this.depthQuarterTarget.texture);
        this.blurFirstPassRenderable = getSsaoBlurRenderable(webgl, this.ssaoDepthTexture, 'horizontal');
        this.blurSecondPassRenderable = getSsaoBlurRenderable(webgl, this.depthBlurProxyTexture, 'vertical');
    }

    setSize(width: number, height: number) {
        const [w, h] = this.texSize;
        const ssaoScale = this.calcSsaoScale(1);
        if (width !== w || height !== h || this.ssaoScale !== ssaoScale) {
            this.texSize.splice(0, 2, width, height);

            const sw = Math.floor(width * this.ssaoScale);
            const sh = Math.floor(height * this.ssaoScale);
            this.downsampledDepthTarget.setSize(sw, sh);
            this.ssaoDepthTexture.define(sw, sh);
            this.depthBlurProxyTexture.define(sw, sh);

            const hw = Math.max(1, Math.floor(sw * 0.5));
            const hh = Math.max(1, Math.floor(sh * 0.5));
            this.depthHalfTarget.setSize(hw, hh);

            const qw = Math.max(1, Math.floor(sw * 0.25));
            const qh = Math.max(1, Math.floor(sh * 0.25));
            this.depthQuarterTarget.setSize(qw, qh);

            ValueCell.update(this.downsampleDepthRenderable.values.uTexSize, Vec2.set(this.downsampleDepthRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.depthHalfRenderable.values.uTexSize, Vec2.set(this.depthHalfRenderable.values.uTexSize.ref.value, hw, hh));
            ValueCell.update(this.depthQuarterRenderable.values.uTexSize, Vec2.set(this.depthQuarterRenderable.values.uTexSize.ref.value, qw, qh));
            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.blurFirstPassRenderable.values.uTexSize, Vec2.set(this.blurFirstPassRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.blurSecondPassRenderable.values.uTexSize, Vec2.set(this.blurSecondPassRenderable.values.uTexSize.ref.value, sw, sh));

            const depthTexture = this.getDepthTexture();
            ValueCell.update(this.depthHalfRenderable.values.tColor, depthTexture);
            ValueCell.update(this.renderable.values.tDepth, depthTexture);

            this.depthHalfRenderable.update();
            this.renderable.update();
        }
    }

    update(camera: ICamera, props: SsaoProps) {
        let needsUpdateSsao = false;
        let needsUpdateSsaoBlur = false;
        let needsUpdateDepthHalf = false;

        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;

        const invProjection = Mat4.identity();
        Mat4.invert(invProjection, camera.projection);

        const [w, h] = this.texSize;
        const v = camera.viewport;

        ValueCell.update(this.renderable.values.uProjection, camera.projection);
        ValueCell.update(this.renderable.values.uInvProjection, invProjection);

        const b = this.renderable.values.uBounds;
        const s = this.ssaoScale;
        Vec4.set(b.ref.value,
            Math.floor(v.x * s) / (w * s),
            Math.floor(v.y * s) / (h * s),
            Math.ceil((v.x + v.width) * s) / (w * s),
            Math.ceil((v.y + v.height) * s) / (h * s)
        );
        ValueCell.update(b, b.ref.value);
        ValueCell.update(this.blurFirstPassRenderable.values.uBounds, b.ref.value);
        ValueCell.update(this.blurSecondPassRenderable.values.uBounds, b.ref.value);

        ValueCell.updateIfChanged(this.blurFirstPassRenderable.values.uNear, camera.near);
        ValueCell.updateIfChanged(this.blurSecondPassRenderable.values.uNear, camera.near);

        ValueCell.updateIfChanged(this.blurFirstPassRenderable.values.uFar, camera.far);
        ValueCell.updateIfChanged(this.blurSecondPassRenderable.values.uFar, camera.far);

        ValueCell.update(this.blurFirstPassRenderable.values.uInvProjection, invProjection);
        ValueCell.update(this.blurSecondPassRenderable.values.uInvProjection, invProjection);

        ValueCell.update(this.blurFirstPassRenderable.values.uBlurDepthBias, props.blurDepthBias);
        ValueCell.update(this.blurSecondPassRenderable.values.uBlurDepthBias, props.blurDepthBias);

        if (this.blurFirstPassRenderable.values.dOrthographic.ref.value !== orthographic) {
            needsUpdateSsaoBlur = true;
            ValueCell.update(this.blurFirstPassRenderable.values.dOrthographic, orthographic);
            ValueCell.update(this.blurSecondPassRenderable.values.dOrthographic, orthographic);
        }

        if (this.nSamples !== props.samples) {
            needsUpdateSsao = true;

            this.nSamples = props.samples;
            ValueCell.update(this.renderable.values.uSamples, getSamples(this.nSamples));
            ValueCell.updateIfChanged(this.renderable.values.dNSamples, this.nSamples);
        }

        const multiScale = props.multiScale.name === 'on';
        if (this.renderable.values.dMultiScale.ref.value !== multiScale) {
            needsUpdateSsao = true;
            ValueCell.update(this.renderable.values.dMultiScale, multiScale);
        }

        if (props.multiScale.name === 'on') {
            const mp = props.multiScale.params;
            if (!deepEqual(this.levels, mp.levels)) {
                needsUpdateSsao = true;

                this.levels = mp.levels;
                const levels = getLevels(mp.levels);
                ValueCell.updateIfChanged(this.renderable.values.dLevels, levels.count);

                ValueCell.update(this.renderable.values.uLevelRadius, levels.radius);
                ValueCell.update(this.renderable.values.uLevelBias, levels.bias);
            }
            ValueCell.updateIfChanged(this.renderable.values.uNearThreshold, mp.nearThreshold);
            ValueCell.updateIfChanged(this.renderable.values.uFarThreshold, mp.farThreshold);
        } else {
            ValueCell.updateIfChanged(this.renderable.values.uRadius, Math.pow(2, props.radius));
        }
        ValueCell.updateIfChanged(this.renderable.values.uBias, props.bias);

        if (this.blurKernelSize !== props.blurKernelSize) {
            needsUpdateSsaoBlur = true;

            this.blurKernelSize = props.blurKernelSize;
            const kernel = getBlurKernel(this.blurKernelSize);

            ValueCell.update(this.blurFirstPassRenderable.values.uKernel, kernel);
            ValueCell.update(this.blurSecondPassRenderable.values.uKernel, kernel);
            ValueCell.update(this.blurFirstPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
            ValueCell.update(this.blurSecondPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
        }

        const ssaoScale = this.calcSsaoScale(props.resolutionScale);
        if (this.ssaoScale !== ssaoScale) {
            needsUpdateSsao = true;
            needsUpdateDepthHalf = true;

            this.ssaoScale = ssaoScale;

            const sw = Math.floor(w * this.ssaoScale);
            const sh = Math.floor(h * this.ssaoScale);
            this.downsampledDepthTarget.setSize(sw, sh);
            this.ssaoDepthTexture.define(sw, sh);
            this.depthBlurProxyTexture.define(sw, sh);

            const hw = Math.floor(sw * 0.5);
            const hh = Math.floor(sh * 0.5);
            this.depthHalfTarget.setSize(hw, hh);

            const qw = Math.floor(sw * 0.25);
            const qh = Math.floor(sh * 0.25);
            this.depthQuarterTarget.setSize(qw, qh);

            const depthTexture = this.getDepthTexture();
            ValueCell.update(this.depthHalfRenderable.values.tColor, depthTexture);
            ValueCell.update(this.renderable.values.tDepth, depthTexture);

            ValueCell.update(this.renderable.values.tDepthHalf, this.depthHalfTarget.texture);
            ValueCell.update(this.renderable.values.tDepthQuarter, this.depthQuarterTarget.texture);

            ValueCell.update(this.downsampleDepthRenderable.values.uTexSize, Vec2.set(this.downsampleDepthRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.depthHalfRenderable.values.uTexSize, Vec2.set(this.depthHalfRenderable.values.uTexSize.ref.value, hw, hh));
            ValueCell.update(this.depthQuarterRenderable.values.uTexSize, Vec2.set(this.depthQuarterRenderable.values.uTexSize.ref.value, qw, qh));
            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.blurFirstPassRenderable.values.uTexSize, Vec2.set(this.blurFirstPassRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.blurSecondPassRenderable.values.uTexSize, Vec2.set(this.blurSecondPassRenderable.values.uTexSize.ref.value, sw, sh));
        }

        if (needsUpdateSsao) {
            this.renderable.update();
        }

        if (needsUpdateSsaoBlur) {
            this.blurFirstPassRenderable.update();
            this.blurSecondPassRenderable.update();
        }

        if (needsUpdateDepthHalf) {
            this.depthHalfRenderable.update();
        }
    }

    render(camera: ICamera) {
        if (isTimingMode) this.webgl.timer.mark('SsaoPass.render');

        const { state } = this.webgl;
        const { x, y, width, height } = camera.viewport;

        const sx = Math.floor(x * this.ssaoScale);
        const sy = Math.floor(y * this.ssaoScale);
        const sw = Math.ceil(width * this.ssaoScale);
        const sh = Math.ceil(height * this.ssaoScale);

        state.viewport(sx, sy, sw, sh);
        state.scissor(sx, sy, sw, sh);

        if (this.ssaoScale < 1) {
            if (isTimingMode) this.webgl.timer.mark('SsaoPass.downsample');
            this.downsampledDepthTarget.bind();
            this.downsampleDepthRenderable.render();
            if (isTimingMode) this.webgl.timer.markEnd('SsaoPass.downsample');
        }

        if (isTimingMode) this.webgl.timer.mark('SsaoPass.half');
        this.depthHalfTarget.bind();
        this.depthHalfRenderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('SsaoPass.half');

        if (isTimingMode) this.webgl.timer.mark('SsaoPass.quarter');
        this.depthQuarterTarget.bind();
        this.depthQuarterRenderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('SsaoPass.quarter');

        this.framebuffer.bind();
        this.renderable.render();

        this.blurFirstPassFramebuffer.bind();
        this.blurFirstPassRenderable.render();

        this.blurSecondPassFramebuffer.bind();
        this.blurSecondPassRenderable.render();
        if (isTimingMode) this.webgl.timer.markEnd('SsaoPass.render');
    }
}

const SsaoSchema = {
    ...QuadSchema,
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tDepthHalf: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tDepthQuarter: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),

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
    uBlurDepthBias: UniformSpec('f'),

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
        uBlurDepthBias: ValueCell.create(0.5),

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
