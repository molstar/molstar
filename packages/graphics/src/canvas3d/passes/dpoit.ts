/**
 * Copyright (c) 2022-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from https://github.com/tsherif/webgl2examples, The MIT License, Copyright © 2017 Tarek Sherif, Shuai Shao
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
import type { Renderbuffer } from '@molstar/graphics/gl/renderbuffer';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { isDebugMode } from '@molstar/core/util/debug';
import type { PassTargets } from './pass-targets.js';

export const BlendBackDpoitSchema = {
  ...QuadSchema,
  tDpoitBackColor: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  uTexSize: UniformSpec('v2'),
};
export const BlendBackDpoitShaderCode = ShaderCode('blend-back-dpoit');
type BlendBackDpoitRenderable = ComputeRenderable<Values<typeof BlendBackDpoitSchema>>;

function getBlendBackDpoitRenderable(ctx: GraphicsContext, dopitBlendBackTexture: Texture): BlendBackDpoitRenderable {
  const values: Values<typeof BlendBackDpoitSchema> = {
    ...QuadValues,
    tDpoitBackColor: ValueCell.create(dopitBlendBackTexture),
    uTexSize: ValueCell.create(Vec2.create(dopitBlendBackTexture.getWidth(), dopitBlendBackTexture.getHeight())),
  };

  const schema = { ...BlendBackDpoitSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', BlendBackDpoitShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

export const EvaluateDpoitSchema = {
  ...QuadSchema,
  tDpoitFrontColor: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  uTexSize: UniformSpec('v2'),
};
export const EvaluateDpoitShaderCode = ShaderCode('evaluate-dpoit');
type EvaluateDpoitRenderable = ComputeRenderable<Values<typeof EvaluateDpoitSchema>>;

function getEvaluateDpoitRenderable(ctx: GraphicsContext, dpoitFrontColorTexture: Texture): EvaluateDpoitRenderable {
  const values: Values<typeof EvaluateDpoitSchema> = {
    ...QuadValues,
    tDpoitFrontColor: ValueCell.create(dpoitFrontColorTexture),
    uTexSize: ValueCell.create(Vec2.create(dpoitFrontColorTexture.getWidth(), dpoitFrontColorTexture.getHeight())),
  };

  const schema = { ...EvaluateDpoitSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', EvaluateDpoitShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

const DpoitCapabilities = ['floatBlend', 'depthTexture', 'multipleRenderTargets'] as const;

export class DpoitPass {
  private readonly DEPTH_CLEAR_VALUE = -99999.0; // NOTE same constant is set in shaders
  private readonly MAX_DEPTH = 1.0;
  private readonly MIN_DEPTH = 0.0;

  private passCount = 0;
  private writeId: number;
  private readId: number;

  private readonly blendBackRenderable: BlendBackDpoitRenderable;
  private readonly renderable: EvaluateDpoitRenderable;

  private readonly depthFramebuffers: Framebuffer[];
  private readonly colorFramebuffers: Framebuffer[];

  private readonly depthTextures: Texture[];
  private readonly colorFrontTextures: Texture[];
  /**
   * the back layer of the current peel, one for both ping-pong sides: the peel doesn't read the previous one, each is
   * blended back (`renderBlendBack`) before the next peel clears it
   */
  private readonly colorBackTexture: Texture;

  private _supported = false;
  get supported() {
    return this._supported;
  }

  getByteCount() {
    if (!this._supported) return 0;
    return (
      this.depthTextures[0].getByteCount() +
      this.depthTextures[1].getByteCount() +
      this.colorFrontTextures[0].getByteCount() +
      this.colorFrontTextures[1].getByteCount() +
      this.colorBackTexture.getByteCount()
    );
  }

  bind(b: PassBuilder) {
    // initialize
    this.passCount = 0;

    b.begin(this.depthFramebuffers[0], {
      label: 'dpoit clear',
      clear: { color: [this.DEPTH_CLEAR_VALUE, this.DEPTH_CLEAR_VALUE, 0, 0] },
    });

    b.begin(this.depthFramebuffers[1], {
      label: 'dpoit clear',
      clear: { color: [-this.MIN_DEPTH, this.MAX_DEPTH, 0, 0] },
    });

    b.begin(this.colorFramebuffers[0], { label: 'dpoit clear', clear: { color: [0, 0, 0, 0] } });

    b.begin(this.colorFramebuffers[1], { label: 'dpoit clear', clear: { color: [0, 0, 0, 0] } });

    b.begin(this.depthFramebuffers[0], { label: 'dpoit transparent' });

    return {
      depth: this.depthTextures[1],
      frontColor: this.colorFrontTextures[1],
    };
  }

  bindDualDepthPeeling(b: PassBuilder) {
    this.readId = this.passCount % 2;
    this.writeId = 1 - this.readId; // ping-pong: 0 or 1

    this.passCount += 1; // increment for next pass

    b.begin(this.depthFramebuffers[this.writeId], {
      label: 'dpoit clear',
      clear: { color: [this.DEPTH_CLEAR_VALUE, this.DEPTH_CLEAR_VALUE, 0, 0] },
    });

    b.begin(this.colorFramebuffers[this.writeId], {
      label: 'dpoit clear',
      clear: { color: [0, 0, 0, 0] },
    });

    b.begin(this.depthFramebuffers[this.writeId], { label: 'dpoit transparent' });

    return {
      depth: this.depthTextures[this.readId],
      frontColor: this.colorFrontTextures[this.readId],
    };
  }

  renderBlendBack(b: PassBuilder) {
    b.quad(this.blendBackRenderable, RenderState.Alpha);
  }

  render(b: PassBuilder) {
    ValueCell.update(this.renderable.values.tDpoitFrontColor, this.colorFrontTextures[this.writeId]);

    this.renderable.update();
    b.quad(this.renderable, RenderState.Premultiplied);
  }

  setSize(width: number, height: number) {
    const [w, h] = this.renderable.values.uTexSize.ref.value;
    if (width !== w || height !== h) {
      for (let i = 0; i < 2; i++) {
        this.depthTextures[i].define(width, height);
        this.colorFrontTextures[i].define(width, height);
      }
      this.colorBackTexture.define(width, height);
      ValueCell.update(
        this.renderable.values.uTexSize,
        Vec2.set(this.renderable.values.uTexSize.ref.value, width, height),
      );
      ValueCell.update(
        this.blendBackRenderable.values.uTexSize,
        Vec2.set(this.blendBackRenderable.values.uTexSize.ref.value, width, height),
      );
    }
  }

  private _init() {
    for (let i = 0; i < 2; i++) {
      // depth
      this.depthFramebuffers[i].setDrawBuffers(3);

      this.colorFrontTextures[i].attachFramebuffer(this.depthFramebuffers[i], 'color0');
      this.colorBackTexture.attachFramebuffer(this.depthFramebuffers[i], 'color1');
      this.depthTextures[i].attachFramebuffer(this.depthFramebuffers[i], 'color2');
      this.stencil?.attachFramebuffer(this.depthFramebuffers[i]);

      // color
      this.colorFramebuffers[i].setDrawBuffers(2);

      this.colorFrontTextures[i].attachFramebuffer(this.colorFramebuffers[i], 'color0');
      this.colorBackTexture.attachFramebuffer(this.colorFramebuffers[i], 'color1');
    }
  }

  /** Why dpoit is not available, if it is not: programs or capabilities missing */
  static getUnsupportedReason(graphics: GraphicsContext): string | undefined {
    if (
      !graphics.isProgramSupported(EvaluateDpoitShaderCode) ||
      !graphics.isProgramSupported(BlendBackDpoitShaderCode)
    ) {
      return 'programs not supported';
    }
    const missing = DpoitCapabilities.filter((k) => !graphics.capabilities[k]);
    if (missing.length) return `missing ${missing.join(', ')}`;
  }

  static isSupported(graphics: GraphicsContext) {
    const reason = DpoitPass.getUnsupportedReason(graphics);
    if (reason && isDebugMode) console.log(`"dpoit" not supported: ${reason}`);
    return !reason;
  }

  constructor(
    targets: PassTargets,
    width: number,
    height: number,
    private stencil: Renderbuffer | null = null,
  ) {
    const { graphics } = targets;
    if (!DpoitPass.isSupported(graphics)) return;

    const {
      resources,
      capabilities: { halfFloatRenderTarget, mixedAttachmentFormats },
    } = graphics;

    // textures

    if (mixedAttachmentFormats) {
      this.depthTextures = [
        targets.texture('image-float32', 'rg', 'float', 'nearest'),
        targets.texture('image-float32', 'rg', 'float', 'nearest'),
      ];

      this.colorFrontTextures = halfFloatRenderTarget
        ? [
            targets.texture('image-float16', 'rgba', 'fp16', 'nearest'),
            targets.texture('image-float16', 'rgba', 'fp16', 'nearest'),
          ]
        : [
            targets.texture('image-float32', 'rgba', 'float', 'nearest'),
            targets.texture('image-float32', 'rgba', 'float', 'nearest'),
          ];

      this.colorBackTexture = halfFloatRenderTarget
        ? targets.texture('image-float16', 'rgba', 'fp16', 'nearest')
        : targets.texture('image-float32', 'rgba', 'float', 'nearest');
    } else {
      // webgl1 requires consistent bit plane counts

      this.depthTextures = [
        targets.texture('image-float32', 'rgba', 'float', 'nearest'),
        targets.texture('image-float32', 'rgba', 'float', 'nearest'),
      ];

      this.colorFrontTextures = [
        targets.texture('image-float32', 'rgba', 'float', 'nearest'),
        targets.texture('image-float32', 'rgba', 'float', 'nearest'),
      ];

      this.colorBackTexture = targets.texture('image-float32', 'rgba', 'float', 'nearest');
    }

    this.depthTextures[0].define(width, height);
    this.depthTextures[1].define(width, height);

    this.colorFrontTextures[0].define(width, height);
    this.colorFrontTextures[1].define(width, height);

    this.colorBackTexture.define(width, height);

    // framebuffers

    this.depthFramebuffers = [resources.framebuffer(), resources.framebuffer()];
    this.colorFramebuffers = [resources.framebuffer(), resources.framebuffer()];

    // renderables

    this.blendBackRenderable = getBlendBackDpoitRenderable(graphics, this.colorBackTexture);
    this.renderable = getEvaluateDpoitRenderable(graphics, this.colorFrontTextures[0]);

    this._supported = true;
    this._init();
  }
}
