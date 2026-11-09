/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { DefineSpec, TextureSpec, UniformSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { Texture } from '@molstar/graphics/gl/texture';
import type { Framebuffer } from '@molstar/graphics/gl/framebuffer';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { Renderer } from '@molstar/graphics/gl/renderer';
import type { Scene } from '@molstar/graphics/gl/scene';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ICamera } from '../camera.js';
import type { PassBuilder } from './pass-builder.js';
import type { PassTargets } from './pass-targets.js';

export const ClipGhostParams = {
  color: PD.Color(Color(0x8c8c8c)),
  alpha: PD.Numeric(0.15, { min: 0, max: 1, step: 0.01 }, { description: 'Opacity of the ghost.' }),
  edgeAlpha: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, { description: 'Opacity of the ghost’s edges.' }),
  resolution: PD.Numeric(
    0.5,
    { min: 0.25, max: 1, step: 0.25 },
    { description: 'Resolution of the ghost relative to the drawing buffer (a ghost is faint, its cost per pixel).' },
  ),
};
export type ClipGhostProps = PD.Values<typeof ClipGhostParams>;

export const ClipGhostOptionParams = PD.MappedStatic(
  'off',
  { on: PD.Group(ClipGhostParams), off: PD.Group({}) },
  { cycle: true, description: 'What the clip objects remove drawn as a faint ghost: a tint and its edges.' },
);
export type ClipGhostOption = PD.Values<{ clipGhost: typeof ClipGhostOptionParams }>['clipGhost'];

/** Of the ghost's distance in view space: a neighbor farther apart than this is across an edge */
const EdgeDistance = 0.05;

/**
 * What the clip objects remove, drawn as a ghost over the image: the removed parts' depth at a lower resolution (the
 * renderer's ghost pass, `Renderer.renderClipGhost`: the clip test inverted), then a tint and its edges where it is
 * nearer than the scene (`composite`), as the marking pass reads the scene's depth
 */
export class ClipGhostPass {
  /** packed depth of the ghost (`depth` variant), cleared to 0.0 (the far plane) where nothing is removed */
  readonly target: RenderTarget;
  private readonly renderable: ClipGhostRenderable;
  /**
   * the image's color alone (`colorTexture`, the target's of the last composite): the depths the composite reads are
   * not attached, as the image's framebuffer may have its opaque depth
   */
  private readonly framebuffer: Framebuffer;
  private colorTexture: Texture | undefined;
  private width: number;
  private height: number;
  private scale = ClipGhostParams.resolution.defaultValue;

  /**
   * whether `option` is on and an object of `scene` has a ghost: clip objects, not opted out (`clipGhost`); a scene
   * planned against (`Canvas3D`'s precompile, its groups counts only) may have one
   */
  static isEnabled(option: ClipGhostOption, scene: Scene) {
    if (option.name !== 'on') return false;
    const { renderables } = scene.primitives;
    for (let i = 0, il = renderables.length; i < il; ++i) {
      const r = renderables[i];
      if (!r) return true;
      if (!r.state.visible || r.state.disposed) continue;
      const values = r.values as { dClipObjects?: ValueCell<boolean>; clipGhost?: ValueCell<boolean> };
      if (values.dClipObjects?.ref.value && values.clipGhost?.ref.value !== false) return true;
    }
    return false;
  }

  constructor(
    private readonly graphics: GraphicsContext,
    targets: PassTargets,
    width: number,
    height: number,
    depth: Texture,
  ) {
    this.width = width;
    this.height = height;
    const [w, h] = this.getTargetSize();
    this.target = targets.createRenderTarget(w, h, 'depth-stencil', 'uint8', 'nearest');
    this.renderable = getClipGhostRenderable(graphics, this.target.texture, depth);
    this.framebuffer = graphics.resources.framebuffer();
  }

  private getTargetSize() {
    return [Math.max(1, Math.round(this.width * this.scale)), Math.max(1, Math.round(this.height * this.scale))];
  }

  getByteCount() {
    return this.target.getByteCount();
  }

  setSize(width: number, height: number) {
    if (width === this.width && height === this.height) return;
    this.width = width;
    this.height = height;
    const [w, h] = this.getTargetSize();
    this.target.setSize(w, h);
  }

  private setScale(scale: number) {
    if (scale === this.scale) return;
    this.scale = scale;
    const [w, h] = this.getTargetSize();
    this.target.setSize(w, h);
  }

  /** the ghost's depth of `scene` as `camera` sees it, at `props.resolution` */
  render(b: PassBuilder, renderer: Renderer, camera: ICamera, scene: Scene, props: ClipGhostProps) {
    // resized when drawn only: a plan (`Canvas3D.prepare` with other props) leaves the target, which the illumination's
    // later iterations composite, as drawn
    if (b.drawing) this.setScale(props.resolution);
    const s = b.drawing ? this.scale : props.resolution;
    b.begin(this.target, { label: 'clip ghost' });
    const { x, y, width, height } = camera.viewport;
    if (b.drawing) {
      renderer.setViewport(
        Math.round(x * s),
        Math.round(y * s),
        Math.max(1, Math.round(width * s)),
        Math.max(1, Math.round(height * s)),
      );
      renderer.clearDepth(true);
    }
    b.renderer({ name: 'clipGhost' }, scene.primitives, () => renderer.renderClipGhost(scene.primitives, camera, s));
    if (b.drawing) renderer.setViewport(x, y, width, height);
  }

  /**
   * The ghost over `target` (the image, of the drawing buffer's size) in `camera`'s viewport, where it is nearer than
   * `depth` (the scene's opaque depth, packed: `packedDepth`) and, if given, the front layer of `transparentDepth`
   * where less transparent than `transparentThreshold`
   */
  composite(
    b: PassBuilder,
    target: RenderTarget,
    camera: ICamera,
    depth: Texture,
    packedDepth: boolean,
    transparentDepth: Texture | null,
    transparentThreshold: number,
    props: ClipGhostProps,
  ) {
    const { values } = this.renderable;
    let changed = false;
    if (values.tDepth.ref.value !== depth) {
      ValueCell.update(values.tDepth, depth);
      changed = true;
    }
    if (values.dPackedDepth.ref.value !== packedDepth) {
      ValueCell.update(values.dPackedDepth, packedDepth);
      changed = true;
    }
    if (transparentDepth && values.tTransparentDepth.ref.value !== transparentDepth) {
      ValueCell.update(values.tTransparentDepth, transparentDepth);
      changed = true;
    }
    if (changed) this.renderable.update();

    if (this.colorTexture !== target.texture) {
      target.texture.attachFramebuffer(this.framebuffer, 'color0');
      this.colorTexture = target.texture;
    }
    b.begin(this.framebuffer, { label: 'clip ghost composite' });
    if (b.drawing) {
      const { x, y, width, height } = camera.viewport;
      this.graphics.setViewport(x, y, width, height);
      this.graphics.setScissor(x, y, width, height);
      ValueCell.update(
        values.uTexSizeInv,
        Vec2.set(values.uTexSizeInv.ref.value, 1 / target.getWidth(), 1 / target.getHeight()),
      );
      ValueCell.update(
        values.uGhostTexSizeInv,
        Vec2.set(values.uGhostTexSizeInv.ref.value, 1 / this.target.getWidth(), 1 / this.target.getHeight()),
      );
      ValueCell.updateIfChanged(values.uHasTransparentDepth, !!transparentDepth);
      ValueCell.updateIfChanged(values.uTransparentThreshold, transparentThreshold);
      ValueCell.updateIfChanged(values.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
      ValueCell.updateIfChanged(values.uNear, camera.near);
      ValueCell.updateIfChanged(values.uFar, camera.far);
      ValueCell.update(values.uColor, Color.toVec3Normalized(values.uColor.ref.value, props.color));
      ValueCell.updateIfChanged(values.uAlpha, props.alpha);
      ValueCell.updateIfChanged(values.uEdgeAlpha, props.edgeAlpha);
    }
    b.quad(this.renderable, RenderState.Premultiplied);
  }

  dispose() {
    this.renderable.dispose();
    this.framebuffer.destroy();
  }
}

//

export const ClipGhostSchema = {
  ...QuadSchema,
  tGhostDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  dPackedDepth: DefineSpec('boolean'),
  tTransparentDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uHasTransparentDepth: UniformSpec('b'),
  uTransparentThreshold: UniformSpec('f'),
  uTexSizeInv: UniformSpec('v2'),
  uGhostTexSizeInv: UniformSpec('v2'),
  uIsOrtho: UniformSpec('f'),
  uNear: UniformSpec('f'),
  uFar: UniformSpec('f'),
  uColor: UniformSpec('v3'),
  uAlpha: UniformSpec('f'),
  uEdgeAlpha: UniformSpec('f'),
  uEdgeDistance: UniformSpec('f'),
};
export const ClipGhostShaderCode = ShaderCode('clip-ghost');
type ClipGhostRenderable = ComputeRenderable<Values<typeof ClipGhostSchema>>;

function getClipGhostRenderable(ctx: GraphicsContext, ghostDepth: Texture, depth: Texture): ClipGhostRenderable {
  const values: Values<typeof ClipGhostSchema> = {
    ...QuadValues,
    tGhostDepth: ValueCell.create(ghostDepth),
    tDepth: ValueCell.create(depth),
    dPackedDepth: ValueCell.create(false),
    // placeholder until a transparent layer is read
    tTransparentDepth: ValueCell.create(ghostDepth),
    uHasTransparentDepth: ValueCell.create(false),
    uTransparentThreshold: ValueCell.create(0.4),
    uTexSizeInv: ValueCell.create(Vec2.create(1, 1)),
    uGhostTexSizeInv: ValueCell.create(Vec2.create(1, 1)),
    uIsOrtho: ValueCell.create(0),
    uNear: ValueCell.create(1),
    uFar: ValueCell.create(10000),
    uColor: ValueCell.create(Vec3()),
    uAlpha: ValueCell.create(0.15),
    uEdgeAlpha: ValueCell.create(0.5),
    uEdgeDistance: ValueCell.create(EdgeDistance),
  };

  const schema = { ...ClipGhostSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', ClipGhostShaderCode, schema, values);
  return createComputeRenderable(renderItem, values);
}
