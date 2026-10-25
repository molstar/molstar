/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { PassBuilder } from './pass-builder.js';
import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { TextureSpec, UniformSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Texture } from '@molstar/graphics/gl/texture';
import { ValueCell } from '@molstar/core/util';
import type { Framebuffer } from '@molstar/graphics/gl/framebuffer';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { isDebugMode } from '@molstar/core/util/debug';
import type { Renderbuffer } from '@molstar/graphics/gl/renderbuffer';
import type { PassTargets } from './pass-targets.js';

export const EvaluateWboitSchema = {
  ...QuadSchema,
  tWboitA: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  tWboitB: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  uTexSize: UniformSpec('v2'),
};
export const EvaluateWboitShaderCode = ShaderCode('evaluate-wboit');
type EvaluateWboitRenderable = ComputeRenderable<Values<typeof EvaluateWboitSchema>>;

function getEvaluateWboitRenderable(
  ctx: GraphicsContext,
  wboitATexture: Texture,
  wboitBTexture: Texture,
): EvaluateWboitRenderable {
  const values: Values<typeof EvaluateWboitSchema> = {
    ...QuadValues,
    tWboitA: ValueCell.create(wboitATexture),
    tWboitB: ValueCell.create(wboitBTexture),
    uTexSize: ValueCell.create(Vec2.create(wboitATexture.getWidth(), wboitATexture.getHeight())),
  };

  const schema = { ...EvaluateWboitSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', EvaluateWboitShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

//

const WboitCapabilities = ['floatBlend', 'depthTexture', 'multipleRenderTargets'] as const;

export class WboitPass {
  private readonly renderable: EvaluateWboitRenderable;

  private readonly framebuffer: Framebuffer;
  private readonly textureA: Texture;
  private readonly textureB: Texture;
  /** the stencil (and depth) of solid interior caps, shared, see the constructor */
  private readonly depthRenderbuffer: Renderbuffer | null;

  private _supported = false;
  get supported() {
    return this._supported;
  }

  getByteCount() {
    if (!this._supported) return 0;
    return this.textureA.getByteCount() + this.textureB.getByteCount();
  }

  bind(b: PassBuilder) {
    b.begin(this.framebuffer, { label: 'wboit transparent', clear: { color: [0, 0, 0, 1] } });
  }

  render(b: PassBuilder) {
    this.renderable.update();
    b.quad(this.renderable, RenderState.Alpha);
  }

  setSize(width: number, height: number) {
    const [w, h] = this.renderable.values.uTexSize.ref.value;
    if (width !== w || height !== h) {
      this.textureA.define(width, height);
      this.textureB.define(width, height);
      ValueCell.update(
        this.renderable.values.uTexSize,
        Vec2.set(this.renderable.values.uTexSize.ref.value, width, height),
      );
    }
  }

  private _init() {
    this.framebuffer.setDrawBuffers(2);

    this.textureA.attachFramebuffer(this.framebuffer, 'color0');
    this.textureB.attachFramebuffer(this.framebuffer, 'color1');

    this.depthRenderbuffer?.attachFramebuffer(this.framebuffer);
  }

  /** Why wboit is not available, if it is not: programs or capabilities missing */
  static getUnsupportedReason(graphics: GraphicsContext): string | undefined {
    if (!graphics.isProgramSupported(EvaluateWboitShaderCode)) return 'programs not supported';
    const missing = WboitCapabilities.filter((k) => !graphics.capabilities[k]);
    if (missing.length) return `missing ${missing.join(', ')}`;
  }

  static isSupported(graphics: GraphicsContext) {
    const reason = WboitPass.getUnsupportedReason(graphics);
    if (reason && isDebugMode) console.log(`"wboit" not supported: ${reason}`);
    return !reason;
  }

  /**
   * `depthStencil`: the depth-stencil buffer of the passes' size the transparent objects are drawn with, only used by
   * solid interior caps (their stencil, cleared per cap) and for the depth the impostors write; shared with a pass
   * that uses it earlier in the frame (the transparent depth target), not of this pass, which owns no depth
   */
  constructor(targets: PassTargets, width: number, height: number, depthStencil: Renderbuffer | null = null) {
    const { graphics } = targets;
    if (!WboitPass.isSupported(graphics)) return;

    const { resources } = graphics;

    const { mixedAttachmentFormats, halfFloatRenderTarget } = graphics.capabilities;
    // the weighted colors and the revealage: half floats, as in the original method, where the attachments may differ
    // in format (WebGL1: not, both float)
    this.textureA =
      mixedAttachmentFormats && halfFloatRenderTarget
        ? targets.texture('image-float16', 'rgba', 'fp16', 'nearest')
        : targets.texture('image-float32', 'rgba', 'float', 'nearest');
    this.textureA.define(width, height);

    // the summed weights (read from `.r`): a red target where the attachments may differ in format
    this.textureB = mixedAttachmentFormats
      ? targets.texture('image-float32', 'r', 'float', 'nearest')
      : targets.texture('image-float32', 'rgba', 'float', 'nearest');
    this.textureB.define(width, height);

    this.depthRenderbuffer = depthStencil;

    this.renderable = getEvaluateWboitRenderable(graphics, this.textureA, this.textureB);
    this.framebuffer = resources.framebuffer();

    this._supported = true;
    this._init();
  }
}
