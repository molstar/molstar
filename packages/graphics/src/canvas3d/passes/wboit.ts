/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

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

export class WboitPass {
  private readonly renderable: EvaluateWboitRenderable;

  private readonly framebuffer: Framebuffer;
  private readonly textureA: Texture;
  private readonly textureB: Texture;
  private readonly depthRenderbuffer: Renderbuffer;

  private _supported = false;
  get supported() {
    return this._supported;
  }

  getByteCount() {
    if (!this._supported) return 0;
    return this.textureA.getByteCount() + this.textureB.getByteCount() + this.depthRenderbuffer.getByteCount();
  }

  bind() {
    this.graphics.beginRenderPass(this.framebuffer, { label: 'wboit transparent', clear: { color: [0, 0, 0, 1] } });
  }

  render() {
    this.graphics.setRenderState(RenderState.Alpha);

    this.renderable.update();
    this.renderable.render();
  }

  setSize(width: number, height: number) {
    const [w, h] = this.renderable.values.uTexSize.ref.value;
    if (width !== w || height !== h) {
      this.textureA.define(width, height);
      this.textureB.define(width, height);
      this.depthRenderbuffer.setSize(width, height);
      ValueCell.update(
        this.renderable.values.uTexSize,
        Vec2.set(this.renderable.values.uTexSize.ref.value, width, height),
      );
    }
  }

  reset() {
    if (this._supported) this._init();
  }

  private _init() {
    this.framebuffer.setDrawBuffers(2);

    this.textureA.attachFramebuffer(this.framebuffer, 'color0');
    this.textureB.attachFramebuffer(this.framebuffer, 'color1');

    this.depthRenderbuffer.attachFramebuffer(this.framebuffer);
  }

  static isSupported(graphics: GraphicsContext) {
    if (!graphics.isProgramSupported(EvaluateWboitShaderCode)) {
      if (isDebugMode) console.log(`Programs required for "wboit" not supported`);
      return false;
    }
    const { floatRenderTarget, depthTexture, multipleRenderTargets } = graphics.capabilities;
    if (!floatRenderTarget || !depthTexture || !multipleRenderTargets) {
      if (isDebugMode) {
        const missing: string[] = [];
        if (!floatRenderTarget) missing.push('floatRenderTarget');
        if (!depthTexture) missing.push('depthTexture');
        if (!multipleRenderTargets) missing.push('multipleRenderTargets');
        console.log(`Missing "${missing.join('", "')}" capabilities required for "wboit"`);
      }
      return false;
    } else {
      return true;
    }
  }

  constructor(
    private graphics: GraphicsContext,
    width: number,
    height: number,
  ) {
    if (!WboitPass.isSupported(graphics)) return;

    const { resources } = graphics;

    this.textureA = resources.texture('image-float32', 'rgba', 'float', 'nearest');
    this.textureA.define(width, height);

    this.textureB = resources.texture('image-float32', 'rgba', 'float', 'nearest');
    this.textureB.define(width, height);

    this.depthRenderbuffer = resources.renderbuffer('depth32f-stencil8', 'depth-stencil', width, height);

    this.renderable = getEvaluateWboitRenderable(graphics, this.textureA, this.textureB);
    this.framebuffer = resources.framebuffer();

    this._supported = true;
    this._init();
  }
}
