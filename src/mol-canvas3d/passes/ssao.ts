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
import { DrawPass } from './draw';
import { ICamera } from '../../mol-canvas3d/camera';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { ssao_frag } from '../../mol-gl/shader/ssao.frag';
import { ssaoDepthMerge_frag } from '../../mol-gl/shader/ssao-depth-merge.frag';
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
    includeOpacity: PD.Boolean(true),
    includeTransparency: PD.Boolean(true),
    separatedTransparency: PD.Boolean(true),
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

    private readonly ssaoFramebuffer: Framebuffer;
    private readonly ssaoBlurFirstPassFramebuffer: Framebuffer;
    private readonly ssaoBlurSecondPassFramebuffer: Framebuffer;

    private readonly downsampledDepthTarget1: RenderTarget;
    private readonly downsampleDepthRenderable1: CopyRenderable;

    private readonly depthHalfTarget1: RenderTarget;
    private readonly depthHalfRenderable1: CopyRenderable;

    private readonly depthQuarterTarget1: RenderTarget;
    private readonly depthQuarterRenderable1: CopyRenderable;

    private readonly downsampledDepthTarget2: RenderTarget;
    private readonly downsampleDepthRenderable2: CopyRenderable;

    private readonly depthHalfTarget2: RenderTarget;
    private readonly depthHalfRenderable2: CopyRenderable;

    private readonly depthQuarterTarget2: RenderTarget;
    private readonly depthQuarterRenderable2: CopyRenderable;

    readonly ssaoDepthTexture: Texture;
    readonly ssaoMergedDepthTarget: RenderTarget;
    // readonly ssaoMergedDepthTexture: Texture;
    readonly ssaoDepthTransparentTexture: Texture;
    private readonly ssaoDepthBlurProxyTexture: Texture;

    private readonly ssaoRenderable: SsaoRenderable;
    private readonly ssaoMergedDephtPassRenderable: SsaoMergeDepthRenderable;
    private readonly ssaoBlurFirstPassRenderable: SsaoBlurRenderable;
    private readonly ssaoBlurSecondPassRenderable: SsaoBlurRenderable;

    private nSamples: number;
    private blurKernelSize: number;
    private texSize: [number, number];
    private includeOpacity: boolean;
    private includeTransparency: boolean;
    private separatedTransparencyPass: boolean;

    private ssaoScale: number;
    private calcSsaoScale(resolutionScale: number) {
        // downscale ssao for high pixel-ratios
        return Math.min(1, 1 / this.webgl.pixelRatio) * resolutionScale;
    }

    private levels: { radius: number, bias: number }[];

    constructor(private readonly webgl: WebGLContext, private readonly drawPass: DrawPass, width: number, height: number) {
        const { textureFloatLinear } = webgl.extensions;
        const { depthTextureOpaque, depthTextureTransparent } = drawPass;

        this.nSamples = 1;
        this.blurKernelSize = 1;
        this.ssaoScale = this.calcSsaoScale(1);
        this.texSize = [width, height];
        this.includeOpacity = true;
        this.includeTransparency = false;
        this.separatedTransparencyPass = false;
        this.levels = [];

        this.ssaoFramebuffer = webgl.resources.framebuffer();
        this.ssaoBlurFirstPassFramebuffer = webgl.resources.framebuffer();
        this.ssaoBlurSecondPassFramebuffer = webgl.resources.framebuffer();

        const sw = Math.floor(width * this.ssaoScale);
        const sh = Math.floor(height * this.ssaoScale);

        const hw = Math.max(1, Math.floor(sw * 0.5));
        const hh = Math.max(1, Math.floor(sh * 0.5));

        const qw = Math.max(1, Math.floor(sw * 0.25));
        const qh = Math.max(1, Math.floor(sh * 0.25));

        const filter = textureFloatLinear ? 'linear' : 'nearest';

        this.downsampledDepthTarget1 = drawPass.packedDepth
            ? webgl.createRenderTarget(sw, sh, false, 'uint8', filter, 'rgba')
            : webgl.createRenderTarget(sw, sh, false, 'float32', filter, webgl.isWebGL2 ? 'alpha' : 'rgba');
        this.downsampleDepthRenderable1 = createCopyRenderable(webgl, depthTextureOpaque);

        const depthTexture1 = this.ssaoScale === 1 ? depthTextureOpaque : this.downsampledDepthTarget1.texture;

        this.depthHalfTarget1 = drawPass.packedDepth
            ? webgl.createRenderTarget(hw, hh, false, 'uint8', filter, 'rgba')
            : webgl.createRenderTarget(hw, hh, false, 'float32', filter, webgl.isWebGL2 ? 'alpha' : 'rgba');
        this.depthHalfRenderable1 = createCopyRenderable(webgl, depthTexture1);

        this.depthQuarterTarget1 = drawPass.packedDepth
            ? webgl.createRenderTarget(qw, qh, false, 'uint8', filter, 'rgba')
            : webgl.createRenderTarget(qw, qh, false, 'float32', filter, webgl.isWebGL2 ? 'alpha' : 'rgba');
        this.depthQuarterRenderable1 = createCopyRenderable(webgl, this.depthHalfTarget1.texture);

        this.downsampledDepthTarget2 = webgl.createRenderTarget(sw, sh, false, 'uint8', filter, 'rgba');
        this.downsampleDepthRenderable2 = createCopyRenderable(webgl, depthTextureTransparent);

        this.ssaoMergedDepthTarget = webgl.createRenderTarget(width, height, false, 'uint8', 'nearest', 'rgba');
        this.ssaoMergedDephtPassRenderable = getSsaoMergeDepthRenderable(webgl, depthTextureOpaque, depthTextureTransparent);

        const depthTexture2 = this.ssaoScale === 1 ? (this.separatedTransparencyPass ? this.ssaoMergedDepthTarget.texture : depthTextureTransparent) : this.downsampledDepthTarget2.texture;

        this.depthHalfTarget2 = webgl.createRenderTarget(hw, hh, false, 'uint8', filter, 'rgba');
        this.depthHalfRenderable2 = createCopyRenderable(webgl, depthTexture2);

        this.depthQuarterTarget2 = webgl.createRenderTarget(qw, qh, false, 'uint8', filter, 'rgba');
        this.depthQuarterRenderable2 = createCopyRenderable(webgl, this.depthHalfTarget2.texture);

        this.ssaoDepthTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', filter);
        this.ssaoDepthTexture.define(sw, sh);
        this.ssaoDepthTexture.attachFramebuffer(this.ssaoFramebuffer, 'color0');

        this.ssaoDepthTransparentTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', filter);
        this.ssaoDepthTransparentTexture.define(sw, sh);

        this.ssaoDepthBlurProxyTexture = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', filter);
        this.ssaoDepthBlurProxyTexture.define(sw, sh);
        this.ssaoDepthBlurProxyTexture.attachFramebuffer(this.ssaoBlurFirstPassFramebuffer, 'color0');

        this.ssaoDepthTexture.attachFramebuffer(this.ssaoBlurSecondPassFramebuffer, 'color0');

        this.ssaoRenderable = getSsaoRenderable(webgl, depthTexture1, this.depthHalfTarget1.texture, this.depthQuarterTarget1.texture, depthTexture2, this.depthHalfTarget2.texture, this.depthQuarterTarget2.texture);
        this.ssaoBlurFirstPassRenderable = getSsaoBlurRenderable(webgl, this.ssaoDepthTexture, 'horizontal');
        this.ssaoBlurSecondPassRenderable = getSsaoBlurRenderable(webgl, this.ssaoDepthBlurProxyTexture, 'vertical');
    }

    setSize(width: number, height: number) {
        const [w, h] = this.texSize;
        const ssaoScale = this.calcSsaoScale(1);
        if (width !== w || height !== h || this.ssaoScale !== ssaoScale) {
            this.texSize.splice(0, 2, width, height);

            ValueCell.update(this.ssaoMergedDephtPassRenderable.values.uTexSize, Vec2.set(this.ssaoMergedDephtPassRenderable.values.uTexSize.ref.value, width, height));
            this.ssaoMergedDepthTarget.setSize(width, height);

            const sw = Math.floor(width * this.ssaoScale);
            const sh = Math.floor(height * this.ssaoScale);
            this.ssaoDepthTexture.define(sw, sh);
            this.ssaoDepthTransparentTexture.define(sw, sh);
            this.ssaoDepthBlurProxyTexture.define(sw, sh);

            const hw = Math.max(1, Math.floor(sw * 0.5));
            const hh = Math.max(1, Math.floor(sh * 0.5));
            const qw = Math.max(1, Math.floor(sw * 0.25));
            const qh = Math.max(1, Math.floor(sh * 0.25));

            this.downsampledDepthTarget1.setSize(sw, sh);
            this.depthHalfTarget1.setSize(hw, hh);
            this.depthQuarterTarget1.setSize(qw, qh);

            ValueCell.update(this.downsampleDepthRenderable1.values.uTexSize, Vec2.set(this.downsampleDepthRenderable1.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.depthHalfRenderable1.values.uTexSize, Vec2.set(this.depthHalfRenderable1.values.uTexSize.ref.value, hw, hh));
            ValueCell.update(this.depthQuarterRenderable1.values.uTexSize, Vec2.set(this.depthQuarterRenderable1.values.uTexSize.ref.value, qw, qh));

            this.downsampledDepthTarget2.setSize(sw, sh);
            this.depthHalfTarget2.setSize(hw, hh);
            this.depthQuarterTarget2.setSize(qw, qh);

            ValueCell.update(this.downsampleDepthRenderable2.values.uTexSize, Vec2.set(this.downsampleDepthRenderable2.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.depthHalfRenderable2.values.uTexSize, Vec2.set(this.depthHalfRenderable2.values.uTexSize.ref.value, hw, hh));
            ValueCell.update(this.depthQuarterRenderable2.values.uTexSize, Vec2.set(this.depthQuarterRenderable2.values.uTexSize.ref.value, qw, qh));

            ValueCell.update(this.ssaoRenderable.values.uTexSize, Vec2.set(this.ssaoRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uTexSize, Vec2.set(this.ssaoBlurFirstPassRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uTexSize, Vec2.set(this.ssaoBlurSecondPassRenderable.values.uTexSize.ref.value, sw, sh));

            const depthTexture1 = this.ssaoScale === 1 ? this.drawPass.depthTextureOpaque : this.downsampledDepthTarget1.texture;
            ValueCell.update(this.depthHalfRenderable1.values.tColor, depthTexture1);
            const depthTexture2 = this.ssaoScale === 1 ? (this.separatedTransparencyPass ? this.drawPass.depthTextureTransparent : this.ssaoMergedDepthTarget.texture) : this.downsampledDepthTarget2.texture;
            ValueCell.update(this.depthHalfRenderable2.values.tColor, depthTexture2);

            ValueCell.update(this.ssaoRenderable.values.tDepth, depthTexture1);
            ValueCell.update(this.ssaoRenderable.values.tDepthTransparent, depthTexture2);

            this.ssaoMergedDephtPassRenderable.update();
            this.depthHalfRenderable1.update();
            this.depthHalfRenderable2.update();
            this.ssaoRenderable.update();
        }
    }

    update(camera: ICamera, props: SsaoProps) {
        let needsUpdateSsao = false;
        let needsUpdateSsaoBlur = false;
        let needsUpdateDepthHalf = false;
        let needsUpdateDownsample2 = false;

        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;

        const invProjection = Mat4.identity();
        Mat4.invert(invProjection, camera.projection);

        const [w, h] = this.texSize;
        const v = camera.viewport;

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

        ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uBlurDepthBias, props.blurDepthBias);
        ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uBlurDepthBias, props.blurDepthBias);

        if (this.ssaoBlurFirstPassRenderable.values.dOrthographic.ref.value !== orthographic) {
            needsUpdateSsaoBlur = true;
            ValueCell.update(this.ssaoBlurFirstPassRenderable.values.dOrthographic, orthographic);
            ValueCell.update(this.ssaoBlurSecondPassRenderable.values.dOrthographic, orthographic);
        }

        const includeOpacity = props.includeOpacity;
        if (this.includeOpacity !== includeOpacity) {
            needsUpdateSsao = true;

            this.includeOpacity = includeOpacity;
            ValueCell.update(this.ssaoRenderable.values.dIncludeOpacity, includeOpacity);
        }

        const includeTransparency = props.includeTransparency;
        if (this.includeTransparency !== includeTransparency) {
            needsUpdateSsao = true;

            this.includeTransparency = includeTransparency;
            ValueCell.update(this.ssaoRenderable.values.dIncludeTransparency, includeTransparency);
        }

        const separatedTransparency = props.separatedTransparency;
        if (this.separatedTransparencyPass !== separatedTransparency) {
            needsUpdateSsao = true;
            needsUpdateDepthHalf = true;
            needsUpdateDownsample2 = true;

            this.separatedTransparencyPass = separatedTransparency;
            ValueCell.update(this.ssaoRenderable.values.dSeparatedTransparency, separatedTransparency);
            ValueCell.update(this.ssaoRenderable.values.tDepthTransparent, this.ssaoScale === 1 ? (separatedTransparency ? this.drawPass.depthTextureTransparent : this.ssaoMergedDepthTarget.texture) : this.downsampledDepthTarget2.texture);
            ValueCell.update(this.downsampleDepthRenderable2.values.tColor, separatedTransparency ? this.drawPass.depthTextureTransparent : this.ssaoMergedDepthTarget.texture);
            ValueCell.update(this.depthHalfRenderable2.values.tColor, separatedTransparency ? this.drawPass.depthTextureTransparent : this.ssaoMergedDepthTarget.texture);
            if (!separatedTransparency) {
                needsUpdateSsaoBlur = true;

                this.ssaoDepthTransparentTexture.attachFramebuffer(this.ssaoFramebuffer, 'color0');
                this.ssaoDepthTransparentTexture.attachFramebuffer(this.ssaoBlurSecondPassFramebuffer, 'color0');
                ValueCell.update(this.ssaoBlurFirstPassRenderable.values.tSsaoDepth, this.ssaoDepthTransparentTexture);
            }
        }

        if (this.nSamples !== props.samples) {
            needsUpdateSsao = true;

            this.nSamples = props.samples;
            ValueCell.update(this.ssaoRenderable.values.uSamples, getSamples(this.nSamples));
            ValueCell.updateIfChanged(this.ssaoRenderable.values.dNSamples, this.nSamples);
        }

        const multiScale = props.multiScale.name === 'on';
        if (this.ssaoRenderable.values.dMultiScale.ref.value !== multiScale) {
            needsUpdateSsao = true;
            ValueCell.update(this.ssaoRenderable.values.dMultiScale, multiScale);
        }

        if (props.multiScale.name === 'on') {
            const mp = props.multiScale.params;
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
            ValueCell.updateIfChanged(this.ssaoRenderable.values.uRadius, Math.pow(2, props.radius));
        }
        ValueCell.updateIfChanged(this.ssaoRenderable.values.uBias, props.bias);

        if (this.blurKernelSize !== props.blurKernelSize) {
            needsUpdateSsaoBlur = true;

            this.blurKernelSize = props.blurKernelSize;
            const kernel = getBlurKernel(this.blurKernelSize);

            ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uKernel, kernel);
            ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uKernel, kernel);
            ValueCell.update(this.ssaoBlurFirstPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
            ValueCell.update(this.ssaoBlurSecondPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
        }

        const ssaoScale = this.calcSsaoScale(props.resolutionScale);
        if (this.ssaoScale !== ssaoScale) {
            needsUpdateSsao = true;
            needsUpdateDepthHalf = true;

            this.ssaoScale = ssaoScale;

            const sw = Math.floor(w * this.ssaoScale);
            const sh = Math.floor(h * this.ssaoScale);
            this.ssaoDepthTexture.define(sw, sh);
            this.ssaoDepthTransparentTexture.define(sw, sh);
            this.ssaoDepthBlurProxyTexture.define(sw, sh);

            const hw = Math.floor(sw * 0.5);
            const hh = Math.floor(sh * 0.5);
            const qw = Math.floor(sw * 0.25);
            const qh = Math.floor(sh * 0.25);

            this.downsampledDepthTarget1.setSize(sw, sh);
            this.depthHalfTarget1.setSize(hw, hh);
            this.depthQuarterTarget1.setSize(qw, qh);

            const depthTexture1 = this.ssaoScale === 1 ? this.drawPass.depthTextureOpaque : this.downsampledDepthTarget1.texture;
            ValueCell.update(this.depthHalfRenderable1.values.tColor, depthTexture1);
            ValueCell.update(this.ssaoRenderable.values.tDepth, depthTexture1);

            ValueCell.update(this.ssaoRenderable.values.tDepthHalf, this.depthHalfTarget1.texture);
            ValueCell.update(this.ssaoRenderable.values.tDepthQuarter, this.depthQuarterTarget1.texture);

            ValueCell.update(this.downsampleDepthRenderable1.values.uTexSize, Vec2.set(this.downsampleDepthRenderable1.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.depthHalfRenderable1.values.uTexSize, Vec2.set(this.depthHalfRenderable1.values.uTexSize.ref.value, hw, hh));
            ValueCell.update(this.depthQuarterRenderable1.values.uTexSize, Vec2.set(this.depthQuarterRenderable1.values.uTexSize.ref.value, qw, qh));

            this.downsampledDepthTarget2.setSize(sw, sh);
            this.depthHalfTarget2.setSize(hw, hh);
            this.depthQuarterTarget2.setSize(qw, qh);

            const depthTexture2 = this.ssaoScale === 1 ? (this.separatedTransparencyPass ? this.drawPass.depthTextureTransparent : this.ssaoMergedDepthTarget.texture) : this.downsampledDepthTarget2.texture;
            ValueCell.update(this.depthHalfRenderable2.values.tColor, depthTexture2);
            ValueCell.update(this.ssaoRenderable.values.tDepthTransparent, depthTexture2);

            ValueCell.update(this.ssaoRenderable.values.tDepthHalfTransparent, this.depthHalfTarget2.texture);
            ValueCell.update(this.ssaoRenderable.values.tDepthQuarterTransparent, this.depthQuarterTarget2.texture);

            ValueCell.update(this.downsampleDepthRenderable2.values.uTexSize, Vec2.set(this.downsampleDepthRenderable2.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.depthHalfRenderable2.values.uTexSize, Vec2.set(this.depthHalfRenderable2.values.uTexSize.ref.value, hw, hh));
            ValueCell.update(this.depthQuarterRenderable2.values.uTexSize, Vec2.set(this.depthQuarterRenderable2.values.uTexSize.ref.value, qw, qh));

            ValueCell.update(this.ssaoRenderable.values.uTexSize, Vec2.set(this.ssaoRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.ssaoBlurFirstPassRenderable.values.uTexSize, Vec2.set(this.ssaoBlurFirstPassRenderable.values.uTexSize.ref.value, sw, sh));
            ValueCell.update(this.ssaoBlurSecondPassRenderable.values.uTexSize, Vec2.set(this.ssaoBlurSecondPassRenderable.values.uTexSize.ref.value, sw, sh));
        }

        if (needsUpdateSsao) {
            this.ssaoRenderable.update();
        }

        if (needsUpdateDownsample2) {
            this.downsampleDepthRenderable2.update();
        }

        if (needsUpdateSsaoBlur) {
            this.ssaoBlurFirstPassRenderable.update();
            this.ssaoBlurSecondPassRenderable.update();
        }

        if (needsUpdateDepthHalf) {
            this.depthHalfRenderable1.update();
            this.depthHalfRenderable2.update();
        }
    }

    render(camera: ICamera) {
        if (isTimingMode) this.webgl.timer.mark('SSAO.render');

        const { state } = this.webgl;
        const { x, y, width, height } = camera.viewport;

        const includeTransparency = this.includeTransparency, includeOpacity = this.includeOpacity, separatedTransparency = this.separatedTransparencyPass;

        if (!separatedTransparency && (includeOpacity || includeTransparency)) {
            state.viewport(x, y, width, height);
            state.scissor(x, y, width, height);

            if (isTimingMode) this.webgl.timer.mark('SSAO.mergedepth');
            this.ssaoMergedDepthTarget.bind();
            this.ssaoMergedDephtPassRenderable.render();
            if (isTimingMode) this.webgl.timer.markEnd('SSAO.mergedepth');
        }

        const sx = Math.floor(x * this.ssaoScale);
        const sy = Math.floor(y * this.ssaoScale);
        const sw = Math.ceil(width * this.ssaoScale);
        const sh = Math.ceil(height * this.ssaoScale);

        state.viewport(sx, sy, sw, sh);
        state.scissor(sx, sy, sw, sh);

        if (this.ssaoScale < 1) {
            if (isTimingMode) this.webgl.timer.mark('SSAO.downsample');
            if (separatedTransparency && includeOpacity) {
                this.downsampledDepthTarget1.bind();
                this.downsampleDepthRenderable1.render();
            }
            if (includeTransparency || !separatedTransparency && includeOpacity) {
                this.downsampledDepthTarget2.bind();
                this.downsampleDepthRenderable2.render();
            }
            if (isTimingMode) this.webgl.timer.markEnd('SSAO.downsample');
        }

        if (isTimingMode) this.webgl.timer.mark('SSAO.half');
        if (separatedTransparency && includeOpacity) {
            this.depthHalfTarget1.bind();
            this.depthHalfRenderable1.render();
        }
        if (includeTransparency || !separatedTransparency && includeOpacity) {
            this.depthHalfTarget2.bind();
            this.depthHalfRenderable2.render();
        }
        if (isTimingMode) this.webgl.timer.markEnd('SSAO.half');

        if (isTimingMode) this.webgl.timer.mark('SSAO.quarter');
        if (separatedTransparency && includeOpacity) {
            this.depthQuarterTarget1.bind();
            this.depthQuarterRenderable1.render();
        }
        if (includeTransparency || !separatedTransparency && includeOpacity) {
            this.depthQuarterTarget2.bind();
            this.depthQuarterRenderable2.render();
        }
        if (isTimingMode) this.webgl.timer.markEnd('SSAO.quarter');
        if (separatedTransparency) {
            if (includeOpacity) {
                this.ssaoDepthTexture.attachFramebuffer(this.ssaoFramebuffer, 'color0');
                ValueCell.update(this.ssaoRenderable.values.uTransparencyFlag, 0);
                this.ssaoFramebuffer.bind();
                this.ssaoRenderable.render();

                ValueCell.update(this.ssaoBlurFirstPassRenderable.values.tSsaoDepth, this.ssaoDepthTexture);
                this.ssaoBlurFirstPassRenderable.update();
                this.ssaoBlurFirstPassFramebuffer.bind();
                this.ssaoBlurFirstPassRenderable.render();

                this.ssaoDepthTexture.attachFramebuffer(this.ssaoBlurSecondPassFramebuffer, 'color0');
                this.ssaoBlurSecondPassFramebuffer.bind();
                this.ssaoBlurSecondPassRenderable.render();
            }
            if (includeTransparency) {
                this.ssaoDepthTransparentTexture.attachFramebuffer(this.ssaoFramebuffer, 'color0');
                ValueCell.update(this.ssaoRenderable.values.uTransparencyFlag, 1);
                this.ssaoFramebuffer.bind();
                this.ssaoRenderable.render();

                ValueCell.update(this.ssaoBlurFirstPassRenderable.values.tSsaoDepth, this.ssaoDepthTransparentTexture);
                this.ssaoBlurFirstPassRenderable.update();
                this.ssaoBlurFirstPassFramebuffer.bind();
                this.ssaoBlurFirstPassRenderable.render();

                this.ssaoDepthTransparentTexture.attachFramebuffer(this.ssaoBlurSecondPassFramebuffer, 'color0');
                this.ssaoBlurSecondPassFramebuffer.bind();
                this.ssaoBlurSecondPassRenderable.render();
            }
        } else {
            this.ssaoFramebuffer.bind();
            this.ssaoRenderable.render();

            this.ssaoBlurFirstPassFramebuffer.bind();
            this.ssaoBlurFirstPassRenderable.render();

            this.ssaoBlurSecondPassFramebuffer.bind();
            this.ssaoBlurSecondPassRenderable.render();
        }
        if (isTimingMode) this.webgl.timer.markEnd('SSAO.render');
    }
}

const SsaoSchema = {
    ...QuadSchema,
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tDepthHalf: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tDepthQuarter: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),

    uTransparencyFlag: UniformSpec('i'),
    dIncludeTransparency: DefineSpec('boolean'),
    dSeparatedTransparency: DefineSpec('boolean'),
    dIncludeOpacity: DefineSpec('boolean'),
    tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tDepthHalfTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tDepthQuarterTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),

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

function getSsaoRenderable(ctx: WebGLContext, depthTexture: Texture, depthHalfTexture: Texture, depthQuarterTexture: Texture, transparentDepthTexture: Texture, transparentDepthHalfTexture: Texture, transparentDepthQuarterTexture: Texture): SsaoRenderable {
    const values: Values<typeof SsaoSchema> = {
        ...QuadValues,
        tDepth: ValueCell.create(depthTexture),
        tDepthHalf: ValueCell.create(depthHalfTexture),
        tDepthQuarter: ValueCell.create(depthQuarterTexture),

        dIncludeTransparency: ValueCell.create(false),
        dSeparatedTransparency: ValueCell.create(false),
        uTransparencyFlag: ValueCell.create(0),
        dIncludeOpacity: ValueCell.create(true),
        tDepthTransparent: ValueCell.create(transparentDepthTexture),
        tDepthHalfTransparent: ValueCell.create(transparentDepthHalfTexture),
        tDepthQuarterTransparent: ValueCell.create(transparentDepthQuarterTexture),

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

const SsaoMergeDephtSchema = {
    ...QuadSchema,
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2')
};

type SsaoMergeDepthRenderable = ComputeRenderable<Values<typeof SsaoMergeDephtSchema>>

function getSsaoMergeDepthRenderable(ctx: WebGLContext, ssaoDepthTexture: Texture, ssaoDepthTransparentTexture: Texture): SsaoMergeDepthRenderable {
    const values: Values<typeof SsaoMergeDephtSchema> = {
        ...QuadValues,
        tDepth: ValueCell.create(ssaoDepthTexture),
        tDepthTransparent: ValueCell.create(ssaoDepthTransparentTexture),
        uTexSize: ValueCell.create(Vec2.create(ssaoDepthTexture.getWidth(), ssaoDepthTexture.getHeight())),
    };
    const schema = { ...SsaoMergeDephtSchema };
    const shaderCode = ShaderCode('ssao_depth_merge', quad_vert, ssaoDepthMerge_frag);
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
