/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { TextureSpec, Values, UniformSpec, DefineSpec } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Texture } from '../../mol-gl/webgl/texture';
import { ValueCell } from '../../mol-util';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { createComputeRenderable, ComputeRenderable } from '../../mol-gl/renderable';
import { Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { DrawPass } from './draw';
import { ICamera } from '../../mol-canvas3d/camera';
import { Scene } from '../../mol-gl/scene';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { postprocessing_frag } from '../../mol-gl/shader/postprocessing.frag';
import { Color } from '../../mol-util/color';
import { FxaaParams, FxaaPass } from './fxaa';
import { SmaaParams, SmaaPass } from './smaa';
import { isTimingMode } from '../../mol-util/debug';
import { BackgroundParams, BackgroundPass } from './background';
import { AssetManager } from '../../mol-util/assets';
import { Light } from '../../mol-gl/renderer';
import { CasParams, CasPass } from './cas';
import { DofPass, DofParams } from './dof';
import { BloomParams } from './bloom';
import { OutlinePass, OutlineProps, OutlineParams } from './outline';
import { ShadowPass, ShadowProps, ShadowParams } from './shadow';
import { SsaoPass, SsaoProps, SsaoParams } from './ssao';


const PostprocessingSchema = {
    ...QuadSchema,
    tSsaoDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tSsaoDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tTransparentColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    dBlendTransparency: DefineSpec('boolean'),
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
    dOcclusionSingleDepth: DefineSpec('boolean'),
    dOcclusionIncludeOpacity: DefineSpec('boolean'),
    dOcclusionIncludeTransparency: DefineSpec('boolean'),
    uOcclusionOffset: UniformSpec('v2'),

    dShadowEnable: DefineSpec('boolean'),

    dOutlineEnable: DefineSpec('boolean'),
    dOutlineScale: DefineSpec('number'),
    dTransparentOutline: DefineSpec('boolean'),
};
type PostprocessingRenderable = ComputeRenderable<Values<typeof PostprocessingSchema>>

function getPostprocessingRenderable(ctx: WebGLContext, colorTexture: Texture, transparentColorTexture: Texture, depthTextureOpaque: Texture, depthTextureTransparent: Texture, shadowsTexture: Texture, outlinesTexture: Texture, ssaoDepthTexture: Texture, ssaoDepthTransparentTexture: Texture, transparentOutline: boolean): PostprocessingRenderable {
    const values: Values<typeof PostprocessingSchema> = {
        ...QuadValues,
        tSsaoDepth: ValueCell.create(ssaoDepthTexture),
        tSsaoDepthTransparent: ValueCell.create(ssaoDepthTransparentTexture),
        tColor: ValueCell.create(colorTexture),
        tTransparentColor: ValueCell.create(transparentColorTexture),
        dBlendTransparency: ValueCell.create(true),
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
        dOcclusionSingleDepth: ValueCell.create(false),
        dOcclusionIncludeOpacity: ValueCell.create(true),
        dOcclusionIncludeTransparency: ValueCell.create(false),
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
        on: PD.Group(SsaoParams),
        off: PD.Group({})
    }, { cycle: true, description: 'Darken occluded crevices with the ambient occlusion effect' }),
    shadow: PD.MappedStatic('off', {
        on: PD.Group(ShadowParams),
        off: PD.Group({})
    }, { cycle: true, description: 'Simplistic shadows' }),
    outline: PD.MappedStatic('off', {
        on: PD.Group(OutlineParams),
        off: PD.Group({})
    }, { cycle: true, description: 'Draw outline around 3D objects' }),
    dof: PD.MappedStatic('off', {
        on: PD.Group(DofParams),
        off: PD.Group({})
    }, { cycle: true, description: 'DOF' }),
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
    bloom: PD.MappedStatic('on', {
        on: PD.Group(BloomParams),
        off: PD.Group({})
    }, { cycle: true, description: 'Bloom' }),
};

export type PostprocessingProps = PD.Values<typeof PostprocessingParams>

export class PostprocessingPass {
    static isEnabled(props: PostprocessingProps) {
        return SsaoPass.isEnabled(props) || ShadowPass.isEnabled(props) || OutlinePass.isEnabled(props) || props.background.variant.name !== 'off';
    }

    static isTransparentDepthRequired(scene: Scene, props: PostprocessingProps) {
        return DofPass.isEnabled(props) || OutlinePass.isEnabled(props) && PostprocessingPass.isTransparentOutlineEnabled(props) || SsaoPass.isEnabled(props) && PostprocessingPass.isTransparentSsaoEnabled(scene, props);
    }

    static isTransparentOutlineEnabled(props: PostprocessingProps) {
        return OutlinePass.isEnabled(props) && ((props.outline.params as OutlineProps).includeTransparent ?? true);
    }

    static isTransparentSsaoEnabled(scene: Scene, props: PostprocessingProps) {
        return SsaoPass.isEnabled(props) && SsaoPass.isTransparentEnabled(scene, props.occlusion.params as SsaoProps);
    }

    static isSsaoEnabled(props: PostprocessingProps) {
        return SsaoPass.isEnabled(props);
    }

    readonly target: RenderTarget;

    private readonly renderable: PostprocessingRenderable;

    readonly ssao: SsaoPass;
    readonly shadow: ShadowPass;
    readonly outline: OutlinePass;
    readonly background: BackgroundPass;

    constructor(private readonly webgl: WebGLContext, assetManager: AssetManager, readonly drawPass: DrawPass) {
        const { colorTarget, transparentColorTarget, depthTextureOpaque, depthTextureTransparent, packedDepth } = drawPass;
        const width = colorTarget.getWidth();
        const height = colorTarget.getHeight();

        // needs to be linear for anti-aliasing pass
        this.target = webgl.createRenderTarget(width, height, false, 'uint8', 'linear');

        this.ssao = new SsaoPass(webgl, width, height, packedDepth, depthTextureOpaque, depthTextureTransparent);
        this.shadow = new ShadowPass(webgl, width, height, depthTextureOpaque);
        this.outline = new OutlinePass(webgl, width, height, depthTextureTransparent, depthTextureOpaque);

        this.renderable = getPostprocessingRenderable(webgl, colorTarget.texture, transparentColorTarget.texture, depthTextureOpaque, depthTextureTransparent, this.shadow.target.texture, this.outline.target.texture, this.ssao.ssaoDepthTexture, this.ssao.ssaoDepthTransparentTexture, true);

        this.background = new BackgroundPass(webgl, assetManager, width, height);
    }

    setSize(width: number, height: number) {
        const [w, h] = this.renderable.values.uTexSize.ref.value;

        if (width !== w || height !== h) {
            this.target.setSize(width, height);
            ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, width, height));
        }

        this.ssao.setSize(width, height);
        this.shadow.setSize(width, height);
        this.outline.setSize(width, height);
        this.background.setSize(width, height);
    }

    updateState(camera: ICamera, scene: Scene, transparentBackground: boolean, backgroundColor: Color, props: PostprocessingProps, light: Light, ambientColor: Vec3) {
        let needsUpdateMain = false;

        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
        const outlinesEnabled = OutlinePass.isEnabled(props);
        const shadowsEnabled = ShadowPass.isEnabled(props);
        const occlusionEnabled = SsaoPass.isEnabled(props);

        if (occlusionEnabled) {
            const params = props.occlusion.params as SsaoProps;
            this.ssao.update(camera, scene, params);
            const includeTransparency = SsaoPass.isTransparentEnabled(scene, params);
            if (this.renderable.values.dOcclusionIncludeTransparency.ref.value !== includeTransparency) {
                needsUpdateMain = true;
                ValueCell.update(this.renderable.values.dOcclusionIncludeTransparency, includeTransparency);
            }
            ValueCell.update(this.renderable.values.uOcclusionColor, Color.toVec3Normalized(this.renderable.values.uOcclusionColor.ref.value, params.color));
        }

        if (shadowsEnabled) {
            this.shadow.update(camera, light, ambientColor, props.shadow.params as ShadowProps);
        }

        if (outlinesEnabled) {
            const outlineProps = props.outline.params as OutlineProps;
            const { transparentOutline, outlineScale } = this.outline.update(camera, outlineProps, this.drawPass.depthTextureTransparent, this.drawPass.depthTextureOpaque);

            ValueCell.update(this.renderable.values.uOutlineColor, Color.toVec3Normalized(this.renderable.values.uOutlineColor.ref.value, outlineProps.color));

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

        const blendTransparency = scene.opacityAverage < 1;
        if (this.renderable.values.dBlendTransparency.ref.value !== blendTransparency) {
            needsUpdateMain = true;
            ValueCell.update(this.renderable.values.dBlendTransparency, blendTransparency);
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

    render(camera: ICamera, scene: Scene, toDrawingBuffer: boolean, transparentBackground: boolean, backgroundColor: Color, props: PostprocessingProps, light: Light, ambientColor: Vec3) {
        if (isTimingMode) this.webgl.timer.mark('PostprocessingPass.render');
        this.updateState(camera, scene, transparentBackground, backgroundColor, props, light, ambientColor);

        const { state } = this.webgl;
        const { x, y, width, height } = camera.viewport;

        // don't render occlusion if offset is given,
        // which will reuse the existing occlusion
        if (props.occlusion.name === 'on' && this.occlusionOffset[0] === 0 && this.occlusionOffset[1] === 0) {
            this.ssao.render(camera);
        }

        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        if (props.outline.name === 'on') {
            this.outline.render();
        }

        if (props.shadow.name === 'on') {
            this.shadow.render();
        }

        if (toDrawingBuffer) {
            this.webgl.unbindFramebuffer();
        } else {
            this.target.bind();
        }

        this.background.update(camera, props.background);
        this.background.clear(props.background, this.transparentBackground, backgroundColor);
        this.background.render(props.background);

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

    constructor(webgl: WebGLContext, width: number, height: number) {
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

    private _renderFxaa(camera: ICamera, input: Texture, target: RenderTarget | undefined, props: PostprocessingProps) {
        if (props.antialiasing.name !== 'fxaa') return;

        this.fxaa.update(input, props.antialiasing.params);
        this.fxaa.render(camera.viewport, target);
    }

    private _renderSmaa(camera: ICamera, input: Texture, target: RenderTarget | undefined, props: PostprocessingProps) {
        if (props.antialiasing.name !== 'smaa') return;

        this.smaa.update(input, props.antialiasing.params);
        this.smaa.render(camera.viewport, target);
    }

    private _renderAntialiasing(camera: ICamera, input: Texture, target: RenderTarget | undefined, props: PostprocessingProps) {
        if (props.antialiasing.name === 'fxaa') {
            this._renderFxaa(camera, input, target, props);
        } else if (props.antialiasing.name === 'smaa') {
            this._renderSmaa(camera, input, target, props);
        }
    }

    private _renderCas(camera: ICamera, input: Texture, target: RenderTarget | undefined, props: PostprocessingProps) {
        if (props.sharpening.name !== 'on') return;

        if (props.antialiasing.name !== 'off') input = this.internalTarget.texture;
        this.cas.update(input, props.sharpening.params);
        this.cas.render(camera.viewport, target);
    }

    render(camera: ICamera, input: Texture, toDrawingBuffer: boolean | RenderTarget, props: PostprocessingProps) {
        if (props.antialiasing.name === 'off' && props.sharpening.name === 'off') return;

        if (props.antialiasing.name === 'smaa' && !this.smaa.supported) {
            console.error('SMAA not supported, missing "HTMLImageElement"');
            return;
        }

        const target = toDrawingBuffer === true
            ? undefined : toDrawingBuffer === false
                ? this.target : toDrawingBuffer;
        if (props.sharpening.name === 'off') {
            this._renderAntialiasing(camera, input, target, props);
        } else if (props.antialiasing.name === 'off') {
            this._renderCas(camera, input, target, props);
        } else {
            this._renderAntialiasing(camera, input, this.internalTarget, props);
            this._renderCas(camera, input, target, props);
        }
    }
}
