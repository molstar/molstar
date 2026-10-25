/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type {
  GraphicsContext,
  RenderPassTarget,
  RenderTargetDepthStencil,
  RenderTargetFormat,
  RenderTargetType,
} from '@molstar/graphics/gl/context';
import type { Renderbuffer, RenderbufferAttachment, RenderbufferFormat } from '@molstar/graphics/gl/renderbuffer';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { Texture, TextureFilter, TextureFormat, TextureKind, TextureType } from '@molstar/graphics/gl/texture';
import { createRenderBuilder, type PassBuilder } from './pass-builder.js';

/** Time without changes or renders before targets are released, see `PassTargets.scheduleRelease` */
export const ReleaseDelayMs = 1000;

function textureKey(id: number) {
  return `texture:${id}`;
}

function renderbufferKey(id: number) {
  return `renderbuffer:${id}`;
}

/** The keys of the textures and renderbuffers `target` draws into */
function getTargetKeys(target: RenderPassTarget): string[] {
  if (target === 'drawing-buffer') return [];
  const framebuffer = 'framebuffer' in target ? target.framebuffer : target;
  return framebuffer.getAttachments().map((r) => `${r.kind}:${r.id}`);
}

/**
 * The render targets of a set of passes (`Passes`, an `ImagePass`), created through it. Their storage is allocated
 * when first drawn into (lazy storage); `releaseUnused` frees the storage of those that neither a plan of the passes
 * for the current props (`createTargetCollector`) nor a frame since the last release drew into, e.g. the targets of
 * a transparency mode or an effect turned off. Data textures (images, lookup tables) are not created through it.
 */
export class PassTargets {
  private readonly textures = new Map<string, Texture>();
  private readonly renderbuffers = new Map<string, Renderbuffer>();
  /** the targets render passes began since the last `scheduleRelease` (frames of the current props and scene) */
  private readonly drawn = new Set<string>();
  private releaseTimer: ReturnType<typeof setTimeout> | undefined;

  constructor(readonly graphics: GraphicsContext) {}

  private addTexture(texture: Texture) {
    const key = textureKey(texture.id);
    this.textures.set(key, texture);
    const destroy = texture.destroy;
    texture.destroy = () => {
      this.textures.delete(key);
      destroy();
    };
    return texture;
  }

  private addRenderbuffer(renderbuffer: Renderbuffer) {
    const key = renderbufferKey(renderbuffer.id);
    this.renderbuffers.set(key, renderbuffer);
    const destroy = renderbuffer.destroy;
    renderbuffer.destroy = () => {
      this.renderbuffers.delete(key);
      destroy();
    };
    return renderbuffer;
  }

  /** See `GraphicsContext.createRenderTarget` */
  createRenderTarget(
    width: number,
    height: number,
    depthStencil?: RenderTargetDepthStencil,
    type?: RenderTargetType,
    filter?: TextureFilter,
    format?: RenderTargetFormat,
  ): RenderTarget {
    const target = this.graphics.createRenderTarget(width, height, depthStencil, type, filter, format);
    this.addTexture(target.texture);
    if (target.depthRenderbuffer) this.addRenderbuffer(target.depthRenderbuffer);
    return target;
  }

  /** See `GraphicsResources.texture`, for textures rendered into (not data) */
  texture(kind: TextureKind, format: TextureFormat, type: TextureType, filter: TextureFilter): Texture {
    return this.addTexture(this.graphics.resources.texture(kind, format, type, filter));
  }

  /** See `GraphicsResources.renderbuffer` */
  renderbuffer(
    format: RenderbufferFormat,
    attachment: RenderbufferAttachment,
    width: number,
    height: number,
  ): Renderbuffer {
    return this.addRenderbuffer(this.graphics.resources.renderbuffer(format, attachment, width, height));
  }

  /** A render builder (`createRenderBuilder`) that records the targets it begins, see `releaseUnused` */
  createRenderBuilder(): PassBuilder {
    const builder = createRenderBuilder(this.graphics);
    return {
      ...builder,
      begin: (target, options) => {
        for (const key of getTargetKeys(target)) this.drawn.add(key);
        builder.begin(target, options);
      },
    };
  }

  /**
   * Frees the storage of the targets neither in `used` (the keys `createTargetCollector` records) nor drawn into
   * since the last `scheduleRelease`; they are allocated again when drawn into
   */
  releaseUnused(used: ReadonlySet<string>) {
    for (const [key, texture] of this.textures) {
      if (!used.has(key) && !this.drawn.has(key)) texture.releaseStorage();
    }
    for (const [key, renderbuffer] of this.renderbuffers) {
      if (!used.has(key) && !this.drawn.has(key)) renderbuffer.releaseStorage();
    }
    this.drawn.clear();
  }

  /** Frees the storage of all targets (e.g. an image pass between renders) */
  releaseAll() {
    for (const texture of this.textures.values()) texture.releaseStorage();
    for (const renderbuffer of this.renderbuffers.values()) renderbuffer.releaseStorage();
    this.drawn.clear();
  }

  /**
   * Calls `release` (e.g. `releaseUnused` with a plan's targets) after `delayMs` without another call, so that a
   * sequence of changes or renders releases once, after it
   */
  scheduleRelease(release: () => void, delayMs = ReleaseDelayMs) {
    this.cancelRelease();
    // what frames draw from now on counts, not what they drew before the change
    this.drawn.clear();
    this.releaseTimer = setTimeout(() => {
      this.releaseTimer = undefined;
      release();
    }, delayMs);
  }

  cancelRelease() {
    if (this.releaseTimer === undefined) return;
    clearTimeout(this.releaseTimer);
    this.releaseTimer = undefined;
  }

  /** Bytes of the targets' storage */
  getByteCount() {
    let count = 0;
    for (const texture of this.textures.values()) count += texture.getByteCount();
    for (const renderbuffer of this.renderbuffers.values()) count += renderbuffer.getByteCount();
    return count;
  }
}

/**
 * Records the targets a plan begins (as `PassTargets.releaseUnused` takes them), without drawing or preparing
 * anything
 */
export function createTargetCollector(): PassBuilder & { readonly used: Set<string> } {
  const used = new Set<string>();
  return {
    drawing: false,
    used,
    begin: (target) => {
      for (const key of getTargetKeys(target)) used.add(key);
    },
    renderer: () => {},
    quad: () => {},
  };
}
