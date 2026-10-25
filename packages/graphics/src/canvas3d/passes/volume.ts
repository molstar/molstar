/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { TextureSpec, UniformSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { Renderer } from '@molstar/graphics/gl/renderer';
import type { Scene } from '@molstar/graphics/gl/scene';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import { createNullTexture, type Texture } from '@molstar/graphics/gl/texture';
import type { ICamera } from '../camera.js';
import type { PassBuilder } from './pass-builder.js';
import type { PassTargets } from './pass-targets.js';
import type { TransparencyMode } from './draw.js';
import { WboitPass } from './wboit.js';

export const VolumeCompositeSchema = {
  ...QuadSchema,
  tVolume: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uTexSize: UniformSpec('v2'),
  uVolumeSize: UniformSpec('v2'),
  uNear: UniformSpec('f'),
  uFar: UniformSpec('f'),
  uIsOrtho: UniformSpec('f'),
};
export const VolumeCompositeShaderCode = ShaderCode('volume-composite');
type VolumeCompositeRenderable = ComputeRenderable<Values<typeof VolumeCompositeSchema>>;

function createVolumeCompositeRenderable(graphics: GraphicsContext, volume: Texture): VolumeCompositeRenderable {
  const values: Values<typeof VolumeCompositeSchema> = {
    ...QuadValues,
    tVolume: ValueCell.create(volume),
    tDepth: ValueCell.create(createNullTexture()),
    uTexSize: ValueCell.create(Vec2()),
    uVolumeSize: ValueCell.create(Vec2()),
    uNear: ValueCell.create(1),
    uFar: ValueCell.create(10000),
    uIsOrtho: ValueCell.create(0),
  };
  const schema = { ...VolumeCompositeSchema };
  const renderItem = graphics.createComputeRenderItem('triangles', VolumeCompositeShaderCode, schema, values);
  return createComputeRenderable(renderItem, values);
}

/**
 * Direct volumes at a lower resolution than the target (`scale`), composited onto it with depth-aware upsampling
 * (`volume-composite`): rays per pixel are their cost, and volumes are smooth
 */
export class VolumePass {
  private readonly target: RenderTarget;
  private readonly composite: VolumeCompositeRenderable;
  /** weighted blended transparency at the lower resolution: volumes are drawn as in the `wboit` mode */
  private readonly wboit: WboitPass;
  private width = 0;
  private height = 0;
  private scale = 1;

  constructor(targets: PassTargets) {
    this.target = targets.createRenderTarget(1, 1, 'none', 'uint8', 'nearest', 'rgba');
    this.composite = createVolumeCompositeRenderable(targets.graphics, this.target.texture);
    this.wboit = new WboitPass(targets, 1, 1);
  }

  /** Whether volumes are drawn by this pass: below the full resolution, and with a depth texture to upsample with */
  static isEnabled(scale: number, packedDepth: boolean) {
    return scale < 1 && !packedDepth;
  }

  getByteCount() {
    return this.target.getByteCount() + this.wboit.getByteCount();
  }

  private setSize(width: number, height: number, scale: number) {
    if (width === this.width && height === this.height && scale === this.scale) return;
    this.width = width;
    this.height = height;
    this.scale = scale;
    const w = Math.max(1, Math.floor(width * scale));
    const h = Math.max(1, Math.floor(height * scale));
    this.target.setSize(w, h);
    if (this.wboit.supported) this.wboit.setSize(w, h);
    ValueCell.update(this.composite.values.uTexSize, Vec2.set(this.composite.values.uTexSize.ref.value, width, height));
    ValueCell.update(this.composite.values.uVolumeSize, Vec2.set(this.composite.values.uVolumeSize.ref.value, w, h));
  }

  /**
   * Draws the volumes of `scene` at `scale` of `out`'s size, as `transparency` draws them at the full size (weighted
   * blended in the `wboit` mode, the volume pass otherwise), and composites them onto `out` (premultiplied); the opaque
   * depth (`depthTexture`, as large as `out`) must not be attached to `out`
   */
  render(
    b: PassBuilder,
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    depthTexture: Texture,
    out: RenderTarget,
    scale: number,
    transparency: TransparencyMode,
  ) {
    const width = out.getWidth();
    const height = out.getHeight();
    this.setSize(width, height, scale);
    const { x, y, width: vw, height: vh } = camera.viewport;

    const viewport = () => {
      // the volume shader reads the opaque depth at its fragment's position over the drawing buffer size
      renderer.setDrawingBufferSize(this.target.getWidth(), this.target.getHeight());
      renderer.setViewport(
        Math.floor(x * scale),
        Math.floor(y * scale),
        Math.ceil((x + vw) * scale) - Math.floor(x * scale),
        Math.ceil((y + vh) * scale) - Math.floor(y * scale),
      );
    };

    if (transparency === 'wboit' && this.wboit.supported) {
      this.wboit.bind(b);
      if (b.drawing) viewport();
      b.renderer({ name: 'wboitTransparent' }, scene.volumes, () =>
        renderer.renderWboitTransparent(scene.volumes, camera, depthTexture),
      );
      // evaluated over transparent black: a premultiplied layer
      b.begin(this.target, { label: 'volume', clear: { color: [0, 0, 0, 0] } });
      this.wboit.render(b);
    } else {
      b.begin(this.target, { label: 'volume', clear: { color: [0, 0, 0, 0] } });
      if (b.drawing) viewport();
      b.renderer({ name: 'volume' }, scene.volumes, () => renderer.renderVolume(scene.volumes, camera, depthTexture));
    }

    const { values } = this.composite;
    ValueCell.update(values.tDepth, depthTexture);
    ValueCell.updateIfChanged(values.uNear, camera.near);
    ValueCell.updateIfChanged(values.uFar, camera.far);
    ValueCell.updateIfChanged(values.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
    this.composite.update();

    b.begin(out, { label: 'volume composite' });
    if (b.drawing) {
      renderer.setDrawingBufferSize(width, height);
      renderer.setViewport(x, y, vw, vh);
    }
    b.quad(this.composite, RenderState.Premultiplied);
  }
}
