/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PassBuilder } from './pass-builder.js';
import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { TextureSpec, UniformSpec, DefineSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import { type Texture, createNullTexture } from '@molstar/graphics/gl/texture';
import { Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Viewport } from '../camera/util.js';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { ICamera } from '../camera.js';
import type { Sphere3D } from '@molstar/core/math/geometry';
import type { PostprocessingProps } from './postprocessing.js';
import type { PassTargets } from './pass-targets.js';

export const DofParams = {
  blurSize: PD.Numeric(9, { min: 1, max: 32, step: 1 }),
  blurSpread: PD.Numeric(1.0, { min: 0.0, max: 10.0, step: 0.1 }),
  inFocus: PD.Numeric(
    0.0,
    { min: -5000.0, max: 5000.0, step: 1.0 },
    { description: 'Distance from the scene center that will be in focus' },
  ),
  PPM: PD.Numeric(
    20.0,
    { min: 0.0, max: 5000.0, step: 0.1 },
    { description: 'Size of the area that will be in focus' },
  ),
  center: PD.Select('camera-target', PD.arrayToOptions(['scene-center', 'camera-target'])),
  mode: PD.Select('plane', PD.arrayToOptions(['plane', 'sphere'])),
};

export type DofProps = PD.Values<typeof DofParams>;

export class DofPass {
  private readonly graphics: GraphicsContext;
  static isEnabled(props: PostprocessingProps) {
    return props.enabled && props.dof.name !== 'off';
  }

  readonly target: RenderTarget;
  /** the weighted color sums of the rows (`dDofPass` 'horizontal'), with half-float targets only */
  private readonly blurTarget: RenderTarget | undefined;
  private readonly horizontalRenderable: DofRenderable | undefined;
  /** 'vertical' after the rows' sums, else 'single' */
  private readonly renderable: DofRenderable;
  private readonly renderables: DofRenderable[];

  constructor(targets: PassTargets, width: number, height: number) {
    const graphics = (this.graphics = targets.graphics);
    this.target = targets.createRenderTarget(width, height, 'none');

    const nullTexture = createNullTexture();
    // the blur separates (its weights are per sample): 2 x blurSize samples instead of blurSize^2, where the sums of
    // the rows can be kept in half floats
    if (graphics.capabilities.halfFloatRenderTarget) {
      this.blurTarget = targets.createRenderTarget(width, height, 'none', 'fp16', 'nearest');
      this.horizontalRenderable = getDofRenderable(graphics, nullTexture, nullTexture, nullTexture, 'horizontal');
      this.renderable = getDofRenderable(graphics, nullTexture, nullTexture, nullTexture, 'vertical');
      ValueCell.update(this.renderable.values.tBlur, this.blurTarget.texture);
      this.renderables = [this.horizontalRenderable, this.renderable];
    } else {
      this.renderable = getDofRenderable(graphics, nullTexture, nullTexture, nullTexture, 'single');
      this.renderables = [this.renderable];
    }
  }

  getByteCount() {
    return this.target.getByteCount() + (this.blurTarget?.getByteCount() ?? 0);
  }

  private updateState(viewport: Viewport) {
    const { x, y, width, height } = viewport;
    this.graphics.setRenderState(RenderState.Default);
    this.graphics.setViewport(x, y, width, height);
    this.graphics.setScissor(x, y, width, height);
    this.graphics.clear({ color: [0, 0, 0, 1] });
  }

  setSize(width: number, height: number) {
    const w = this.target.texture.getWidth();
    const h = this.target.texture.getHeight();

    if (width !== w || height !== h) {
      this.target.setSize(width, height);
      this.blurTarget?.setSize(width, height);
      for (const r of this.renderables) {
        ValueCell.update(r.values.uTexSize, Vec2.set(r.values.uTexSize.ref.value, width, height));
      }
    }
  }

  update(
    camera: ICamera,
    input: Texture,
    depthOpaque: Texture,
    depthTransparent: Texture,
    props: DofProps,
    sphere: Sphere3D,
  ) {
    for (const renderable of this.renderables) {
      this.updateRenderable(renderable, camera, input, depthOpaque, depthTransparent, props, sphere);
    }
  }

  private updateRenderable(
    renderable: DofRenderable,
    camera: ICamera,
    input: Texture,
    depthOpaque: Texture,
    depthTransparent: Texture,
    props: DofProps,
    sphere: Sphere3D,
  ) {
    let needsUpdate = false;
    if (renderable.values.tColor.ref.value !== input) {
      ValueCell.update(renderable.values.tColor, input);
      needsUpdate = true;
    }
    if (renderable.values.tDepthOpaque.ref.value !== depthOpaque) {
      ValueCell.update(renderable.values.tDepthOpaque, depthOpaque);
      needsUpdate = true;
    }
    if (renderable.values.tDepthTransparent.ref.value !== depthTransparent) {
      ValueCell.update(renderable.values.tDepthTransparent, depthTransparent);
      needsUpdate = true;
    }
    const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
    const invProjection = renderable.values.uInvProjection.ref.value;
    Mat4.invert(invProjection, camera.projection);

    const [w, h] = renderable.values.uTexSize.ref.value;
    const v = camera.viewport;

    ValueCell.update(renderable.values.uProjection, camera.projection);
    ValueCell.update(renderable.values.uInvProjection, invProjection);
    ValueCell.update(renderable.values.uMode, props.mode === 'sphere' ? 1 : 0);

    Vec4.set(renderable.values.uBounds.ref.value, v.x / w, v.y / h, (v.x + v.width) / w, (v.y + v.height) / h);
    ValueCell.update(renderable.values.uBounds, renderable.values.uBounds.ref.value);

    ValueCell.updateIfChanged(renderable.values.uNear, camera.near);
    ValueCell.updateIfChanged(renderable.values.uFar, camera.far);
    ValueCell.updateIfChanged(renderable.values.dOrthographic, orthographic);

    const blurSize = Math.round(props.blurSize * this.graphics.pixelRatio);
    if (renderable.values.dBlurSize.ref.value !== blurSize) {
      ValueCell.update(renderable.values.dBlurSize, blurSize);
      needsUpdate = true;
    }

    const worldCenter = props.center === 'scene-center' ? sphere.center : camera.state.target;
    const distance = Vec3.distance(camera.state.position, worldCenter);
    const inFocus = distance + props.inFocus;
    ValueCell.updateIfChanged(renderable.values.uInFocus, inFocus * camera.scale);

    // transform center in view space
    const center = renderable.values.uCenter.ref.value;
    Vec3.transformMat4(center, worldCenter, camera.view);
    ValueCell.update(renderable.values.uCenter, center);

    ValueCell.updateIfChanged(renderable.values.uBlurSpread, props.blurSpread);
    ValueCell.updateIfChanged(renderable.values.uPPM, props.PPM * camera.scale);

    if (needsUpdate) {
      renderable.update();
    }
  }

  render(b: PassBuilder, viewport: Viewport, target: undefined | RenderTarget) {
    if (this.blurTarget && this.horizontalRenderable) {
      b.begin(this.blurTarget, { label: 'dof horizontal' });
      if (b.drawing) this.updateState(viewport);
      b.quad(this.horizontalRenderable, RenderState.Default);
    }
    b.begin(target ?? 'drawing-buffer', { label: 'dof' });
    if (b.drawing) this.updateState(viewport);
    b.quad(this.renderable, RenderState.Default);
  }
}

//

export const DofSchema = {
  ...QuadSchema,
  tDepthOpaque: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tBlur: TextureSpec('texture', 'rgba', 'fp16', 'nearest', undefined, { usedIf: ['dDofPass_vertical'] }),
  uTexSize: UniformSpec('v2'),
  dDofPass: DefineSpec('string', ['single', 'horizontal', 'vertical']),

  uProjection: { ...UniformSpec('m4'), depthRange: 'projection' as const },
  uInvProjection: { ...UniformSpec('m4'), depthRange: 'inverse' as const },
  uBounds: { ...UniformSpec('v4'), windowRect: 'bounds' as const },
  uCenter: UniformSpec('v3'),
  uMode: UniformSpec('i'),

  dOrthographic: DefineSpec('number'),
  uNear: UniformSpec('f'),
  uFar: UniformSpec('f'),

  dBlurSize: DefineSpec('number'),
  uBlurSpread: UniformSpec('f'),
  uInFocus: UniformSpec('f'),
  uPPM: UniformSpec('f'),
};

export const DofShaderCode = ShaderCode('dof');
type DofRenderable = ComputeRenderable<Values<typeof DofSchema>>;

function getDofRenderable(
  ctx: GraphicsContext,
  colorTexture: Texture,
  depthTextureOpaque: Texture,
  depthTextureTransparent: Texture,
  pass: 'single' | 'horizontal' | 'vertical',
): DofRenderable {
  const width = colorTexture.getWidth();
  const height = colorTexture.getHeight();

  const values: Values<typeof DofSchema> = {
    ...QuadValues,
    tDepthOpaque: ValueCell.create(depthTextureOpaque),
    tDepthTransparent: ValueCell.create(depthTextureTransparent),
    tColor: ValueCell.create(colorTexture),
    tBlur: ValueCell.create(createNullTexture()),
    uTexSize: ValueCell.create(Vec2.create(width, height)),
    dDofPass: ValueCell.create(pass),

    uProjection: ValueCell.create(Mat4.identity()),
    uInvProjection: ValueCell.create(Mat4.identity()),
    uBounds: ValueCell.create(Vec4()),
    uCenter: ValueCell.create(Vec3()),
    uMode: ValueCell.create(0),

    dOrthographic: ValueCell.create(0),
    uNear: ValueCell.create(1),
    uFar: ValueCell.create(10000),

    dBlurSize: ValueCell.create(5),
    uBlurSpread: ValueCell.create(300.0),
    uInFocus: ValueCell.create(20.0),
    uPPM: ValueCell.create(20.0),
  };

  const schema = { ...DofSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', DofShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}
