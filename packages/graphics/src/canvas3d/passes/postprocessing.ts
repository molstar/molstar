/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { PassBuilder } from './pass-builder.js';
import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { TextureSpec, type Values, UniformSpec, DefineSpec } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { Texture, TextureFilter } from '@molstar/graphics/gl/texture';
import { ValueCell } from '@molstar/core/util';
import { createComputeRenderable, type ComputeRenderable } from '@molstar/graphics/gl/renderable';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { DrawPass } from './draw.js';
import type { ICamera } from '../camera.js';
import type { Scene } from '@molstar/graphics/gl/scene';
import { Color } from '@molstar/core/util/color';
import { FxaaParams, FxaaPass } from './fxaa.js';
import { SmaaParams, SmaaPass } from './smaa.js';
import { isTimingMode } from '@molstar/core/util/debug';
import { BackgroundParams, BackgroundPass } from './background.js';
import type { AssetManager } from '@molstar/core/util/assets';
import type { Light } from '@molstar/graphics/gl/renderer';
import { CasParams, CasPass } from './cas.js';
import { DofPass, DofParams } from './dof.js';
import { BloomParams, BloomPass } from './bloom.js';
import { OutlinePass, type OutlineProps, OutlineParams } from './outline.js';
import { ShadowPass, type ShadowProps, ShadowParams } from './shadow.js';
import { SsaoPass, type SsaoProps, SsaoParams, updateOcclusionUpsample } from './ssao.js';
import type { PassTargets } from './pass-targets.js';

export const PostprocessingSchema = {
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
  tBloom: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  dBloomEnable: DefineSpec('boolean'),
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
  uOcclusionOffset: { ...UniformSpec('v2'), windowRect: 'offset' as const },
  /** the SSAO targets have a lower resolution, see `getSsaoUpsampled` */
  dOcclusionUpsample: DefineSpec('boolean'),
  uOcclusionTexSize: UniformSpec('v2'),

  dShadowEnable: DefineSpec('boolean'),

  dOutlineEnable: DefineSpec('boolean'),
  dOutlineScale: DefineSpec('number'),
  dTransparentOutline: DefineSpec('boolean'),
};
export const PostprocessingShaderCode = ShaderCode('postprocessing');
type PostprocessingRenderable = ComputeRenderable<Values<typeof PostprocessingSchema>>;

function getPostprocessingRenderable(
  ctx: GraphicsContext,
  colorTexture: Texture,
  transparentColorTexture: Texture,
  depthTextureOpaque: Texture,
  depthTextureTransparent: Texture,
  shadowsTexture: Texture,
  outlinesTexture: Texture,
  ssaoDepthTexture: Texture,
  ssaoDepthTransparentTexture: Texture,
  bloomTexture: Texture,
  transparentOutline: boolean,
): PostprocessingRenderable {
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
    tBloom: ValueCell.create(bloomTexture),
    dBloomEnable: ValueCell.create(false),
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
    dOcclusionUpsample: ValueCell.create(false),
    uOcclusionTexSize: ValueCell.create(Vec2()),

    dShadowEnable: ValueCell.create(false),

    dOutlineEnable: ValueCell.create(false),
    dOutlineScale: ValueCell.create(1),
    dTransparentOutline: ValueCell.create(transparentOutline),
  };

  const schema = { ...PostprocessingSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', PostprocessingShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

export const PostprocessingParams = {
  enabled: PD.Boolean(true),
  occlusion: PD.MappedStatic(
    'on',
    {
      on: PD.Group(SsaoParams),
      off: PD.Group({}),
    },
    {
      cycle: true,
      description: 'Darken occluded crevices with the ambient occlusion effect',
      hideIf: (p) => p.enabled === false,
    },
  ),
  shadow: PD.MappedStatic(
    'off',
    {
      on: PD.Group(ShadowParams),
      off: PD.Group({}),
    },
    { cycle: true, description: 'Simplistic shadows', hideIf: (p) => p.enabled === false },
  ),
  outline: PD.MappedStatic(
    'off',
    {
      on: PD.Group(OutlineParams),
      off: PD.Group({}),
    },
    { cycle: true, description: 'Draw outline around 3D objects', hideIf: (p) => p.enabled === false },
  ),
  dof: PD.MappedStatic(
    'off',
    {
      on: PD.Group(DofParams),
      off: PD.Group({}),
    },
    { cycle: true, description: 'DOF', hideIf: (p) => p.enabled === false },
  ),
  antialiasing: PD.MappedStatic(
    'smaa',
    {
      fxaa: PD.Group(FxaaParams),
      smaa: PD.Group(SmaaParams),
      off: PD.Group({}),
    },
    {
      options: [
        ['fxaa', 'FXAA'],
        ['smaa', 'SMAA'],
        ['off', 'Off'],
      ],
      description: 'Smooth pixel edges',
      hideIf: (p) => p.enabled === false,
    },
  ),
  sharpening: PD.MappedStatic(
    'off',
    {
      on: PD.Group(CasParams),
      off: PD.Group({}),
    },
    { cycle: true, description: 'Contrast Adaptive Sharpening', hideIf: (p) => p.enabled === false },
  ),
  background: PD.Group(BackgroundParams, { isFlat: true, hideIf: (p) => p.enabled === false }),
  bloom: PD.MappedStatic(
    'on',
    {
      on: PD.Group(BloomParams),
      off: PD.Group({}),
    },
    { cycle: true, description: 'Bloom', hideIf: (p) => p.enabled === false },
  ),
};

export type PostprocessingProps = PD.Values<typeof PostprocessingParams>;

export class PostprocessingPass {
  /** Whether the pass draws: bloom only counts when it draws in `scene` (`BloomPass.isActive`) */
  static isEnabled(scene: Scene, props: PostprocessingProps) {
    return (
      props.enabled &&
      (SsaoPass.isEnabled(props) ||
        ShadowPass.isEnabled(props) ||
        OutlinePass.isEnabled(props) ||
        props.background.variant.name !== 'off' ||
        BloomPass.isActive(scene, props))
    );
  }

  /**
   * Whether a pass reads the transparent depth. The target is cleared then also for scenes without transparency, so
   * the passes never read the depth of an earlier frame or scene; transparent objects are only drawn into it if any.
   */
  static isTransparentDepthRequired(scene: Scene, props: PostprocessingProps) {
    return (
      props.enabled &&
      (DofPass.isEnabled(props) ||
        (OutlinePass.isEnabled(props) && PostprocessingPass.isTransparentOutlineEnabled(props)) ||
        (SsaoPass.isEnabled(props) && PostprocessingPass.isTransparentSsaoEnabled(scene, props)) ||
        // bloom reads it when it draws, with blended transparency only, see `BloomPass.update`
        (BloomPass.isActive(scene, props) && scene.opacityAverage < 1))
    );
  }

  static isTransparentOutlineEnabled(props: PostprocessingProps) {
    return (
      props.enabled &&
      OutlinePass.isEnabled(props) &&
      ((props.outline.params as OutlineProps).includeTransparent ?? true)
    );
  }

  static isTransparentSsaoEnabled(scene: Scene, props: PostprocessingProps) {
    return (
      props.enabled &&
      SsaoPass.isEnabled(props) &&
      SsaoPass.isTransparentEnabled(scene, props.occlusion.params as SsaoProps)
    );
  }

  static isSsaoEnabled(props: PostprocessingProps) {
    return props.enabled && SsaoPass.isEnabled(props);
  }

  readonly target: RenderTarget;

  private readonly renderable: PostprocessingRenderable;

  readonly ssao: SsaoPass;
  readonly shadow: ShadowPass;
  readonly outline: OutlinePass;
  readonly bloom: BloomPass;
  readonly background: BackgroundPass;

  constructor(
    private readonly graphics: GraphicsContext,
    assetManager: AssetManager,
    readonly drawPass: DrawPass,
  ) {
    const { colorTarget, transparentColorTarget, depthTextureOpaque, depthTextureTransparent, packedDepth } = drawPass;
    const { targets } = drawPass;
    const width = colorTarget.getWidth();
    const height = colorTarget.getHeight();

    // needs to be linear for anti-aliasing pass
    this.target = targets.createRenderTarget(width, height, 'none', 'uint8', 'linear');

    this.ssao = new SsaoPass(targets, width, height, packedDepth, depthTextureOpaque, depthTextureTransparent);
    this.shadow = new ShadowPass(targets, width, height, depthTextureOpaque);
    this.outline = new OutlinePass(targets, width, height, depthTextureTransparent, depthTextureOpaque);
    this.bloom = new BloomPass(targets, width, height);

    this.renderable = getPostprocessingRenderable(
      graphics,
      colorTarget.texture,
      transparentColorTarget.texture,
      depthTextureOpaque,
      depthTextureTransparent,
      this.shadow.target.texture,
      this.outline.target.texture,
      this.ssao.ssaoDepthTexture,
      this.ssao.ssaoDepthTransparentTexture,
      this.bloom.compositeTarget.texture,
      true,
    );

    this.background = new BackgroundPass(graphics, assetManager, width, height);
  }

  getByteCount() {
    return (
      this.target.getByteCount() +
      this.ssao.getByteCount() +
      this.shadow.getByteCount() +
      this.outline.getByteCount() +
      this.bloom.getByteCount()
    );
  }

  setSize(width: number, height: number) {
    const [w, h] = this.renderable.values.uTexSize.ref.value;

    if (width !== w || height !== h) {
      this.target.setSize(width, height);
      ValueCell.update(
        this.renderable.values.uTexSize,
        Vec2.set(this.renderable.values.uTexSize.ref.value, width, height),
      );
    }

    this.ssao.setSize(width, height);
    this.shadow.setSize(width, height);
    this.outline.setSize(width, height);
    this.bloom.setSize(width, height);
    this.background.setSize(width, height);
  }

  updateState(
    camera: ICamera,
    scene: Scene,
    transparentBackground: boolean,
    backgroundColor: Color,
    props: PostprocessingProps,
    light: Light,
    ambientColor: Vec3,
  ) {
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
      const [width, height] = this.renderable.values.uTexSize.ref.value;
      if (updateOcclusionUpsample(this.renderable.values, this.ssao.ssaoDepthTexture, width, height)) {
        needsUpdateMain = true;
      }
      ValueCell.update(
        this.renderable.values.uOcclusionColor,
        Color.toVec3Normalized(this.renderable.values.uOcclusionColor.ref.value, params.color),
      );
    }

    if (shadowsEnabled) {
      this.shadow.update(camera, light, ambientColor, props.shadow.params as ShadowProps);
    }

    if (outlinesEnabled) {
      const outlineProps = props.outline.params as OutlineProps;
      const { transparentOutline, outlineScale } = this.outline.update(
        camera,
        outlineProps,
        this.drawPass.depthTextureTransparent,
        this.drawPass.depthTextureOpaque,
      );

      ValueCell.update(
        this.renderable.values.uOutlineColor,
        Color.toVec3Normalized(this.renderable.values.uOutlineColor.ref.value, outlineProps.color),
      );

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
    ValueCell.update(
      this.renderable.values.uFogColor,
      Color.toVec3Normalized(this.renderable.values.uFogColor.ref.value, backgroundColor),
    );
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
  }

  private occlusionOffset: [x: number, y: number] = [0, 0];
  /** shift in uv (top-left origin) of the occlusion rendered without jitter to the jittered image */
  setOcclusionOffset(x: number, y: number) {
    this.occlusionOffset[0] = x;
    this.occlusionOffset[1] = y;
    ValueCell.update(
      this.renderable.values.uOcclusionOffset,
      Vec2.set(this.renderable.values.uOcclusionOffset.ref.value, x, y),
    );
  }

  private transparentBackground = false;
  setTransparentBackground(value: boolean) {
    this.transparentBackground = value;
  }

  render(
    b: PassBuilder,
    camera: ICamera,
    scene: Scene,
    toDrawingBuffer: boolean,
    transparentBackground: boolean,
    backgroundColor: Color,
    props: PostprocessingProps,
    light: Light,
    ambientColor: Vec3,
    bloomEnable: boolean = false,
  ) {
    const timing = isTimingMode && b.drawing;
    if (timing) this.graphics.timer.group('PostprocessingPass.render');
    this.updateState(camera, scene, transparentBackground, backgroundColor, props, light, ambientColor);

    if (this.renderable.values.dBloomEnable.ref.value !== bloomEnable) {
      ValueCell.update(this.renderable.values.dBloomEnable, bloomEnable);
      this.renderable.update();
    }

    // premultiplied scene composited over the background (cleared color or environment) with
    // premultiplied-over — keeps transparency weight (vs SRC_ALPHA's color*a^2) and lets bloom occlude
    const compositeBackground = this.background.isEnabled(props) || (bloomEnable && !this.transparentBackground);

    const { x, y, width, height } = camera.viewport;

    // don't render occlusion if offset is given,
    // which will reuse the existing occlusion
    if (props.occlusion.name === 'on' && this.occlusionOffset[0] === 0 && this.occlusionOffset[1] === 0) {
      this.ssao.render(b, camera);
    }

    if (b.drawing) {
      this.graphics.setViewport(x, y, width, height);
      this.graphics.setScissor(x, y, width, height);
    }

    if (props.outline.name === 'on') {
      this.outline.render(b);
    }

    if (props.shadow.name === 'on') {
      this.shadow.render(b);
    }

    b.begin(toDrawingBuffer ? 'drawing-buffer' : this.target, { label: 'postprocessing' });

    if (b.drawing) {
      this.background.update(camera, props.background);
      this.background.clear(props.background, this.transparentBackground, backgroundColor);
      this.background.render(b, props.background);
    } else {
      // the props may be hypothetical (`Canvas3D.prepare`): the live background stays
      this.background.prepare(b, props.background);
    }

    b.quad(this.renderable, compositeBackground ? RenderState.Premultiplied : RenderState.Default);
    if (timing) this.graphics.timer.groupEnd('PostprocessingPass.render');
  }
}

export class AntialiasingPass {
  static isEnabled(props: PostprocessingProps) {
    return props.enabled && (props.antialiasing.name !== 'off' || props.sharpening.name !== 'off');
  }

  readonly target: RenderTarget;
  private readonly internalTarget: RenderTarget;

  private readonly fxaa: FxaaPass;
  private readonly smaa: SmaaPass;
  private readonly cas: CasPass;

  constructor(targets: PassTargets, width: number, height: number, filter: TextureFilter = 'nearest') {
    const { graphics } = targets;
    this.target = targets.createRenderTarget(width, height, 'none', 'uint8', filter);
    this.internalTarget = targets.createRenderTarget(width, height, 'none');

    this.fxaa = new FxaaPass(graphics, this.target.texture);
    this.smaa = new SmaaPass(targets, this.target.texture);
    this.cas = new CasPass(graphics, this.target.texture);
  }

  getByteCount() {
    return this.target.getByteCount() + this.internalTarget.getByteCount() + this.smaa.getByteCount();
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

  private _renderFxaa(
    b: PassBuilder,
    camera: ICamera,
    input: Texture,
    target: RenderTarget | undefined,
    props: PostprocessingProps,
  ) {
    if (props.antialiasing.name !== 'fxaa') return;

    this.fxaa.update(input, props.antialiasing.params);
    this.fxaa.render(b, camera.viewport, target);
  }

  private _renderSmaa(
    b: PassBuilder,
    camera: ICamera,
    input: Texture,
    target: RenderTarget | undefined,
    props: PostprocessingProps,
  ) {
    if (props.antialiasing.name !== 'smaa') return;

    this.smaa.update(input, props.antialiasing.params);
    this.smaa.render(b, camera.viewport, target);
  }

  private _renderAntialiasing(
    b: PassBuilder,
    camera: ICamera,
    input: Texture,
    target: RenderTarget | undefined,
    props: PostprocessingProps,
  ) {
    if (props.antialiasing.name === 'fxaa') {
      this._renderFxaa(b, camera, input, target, props);
    } else if (props.antialiasing.name === 'smaa') {
      this._renderSmaa(b, camera, input, target, props);
    }
  }

  private _renderCas(
    b: PassBuilder,
    camera: ICamera,
    input: Texture,
    target: RenderTarget | undefined,
    props: PostprocessingProps,
  ) {
    if (props.sharpening.name !== 'on') return;

    if (props.antialiasing.name !== 'off') input = this.internalTarget.texture;
    this.cas.update(input, props.sharpening.params);
    this.cas.render(b, camera.viewport, target);
  }

  /** Antialiasing without sharpening, returns false if nothing was rendered. */
  renderAntialiasingOnly(
    b: PassBuilder,
    camera: ICamera,
    input: Texture,
    target: RenderTarget,
    props: PostprocessingProps,
  ): boolean {
    if (props.antialiasing.name === 'off') return false;
    if (props.antialiasing.name === 'smaa' && !this.smaa.supported) return false;

    this._renderAntialiasing(b, camera, input, target, props);
    return true;
  }

  render(
    b: PassBuilder,
    camera: ICamera,
    input: Texture,
    toDrawingBuffer: boolean | RenderTarget,
    props: PostprocessingProps,
  ) {
    if (props.antialiasing.name === 'off' && props.sharpening.name === 'off') return;

    if (props.antialiasing.name === 'smaa' && !this.smaa.supported) {
      console.error('SMAA not supported, missing "HTMLImageElement"');
      return;
    }

    const target = toDrawingBuffer === true ? undefined : toDrawingBuffer === false ? this.target : toDrawingBuffer;
    if (props.sharpening.name === 'off') {
      this._renderAntialiasing(b, camera, input, target, props);
    } else if (props.antialiasing.name === 'off') {
      this._renderCas(b, camera, input, target, props);
    } else {
      this._renderAntialiasing(b, camera, input, this.internalTarget, props);
      this._renderCas(b, camera, input, target, props);
    }
  }
}
