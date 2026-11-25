/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsResources } from './context.js';
import type { Renderbuffer } from './renderbuffer.js';
import type { Texture } from './texture.js';

/** A layer of a 3d texture as a color attachment */
export interface RenderPassColorLayer {
  readonly texture: Texture;
  readonly layer: number;
}

export type RenderPassColorAttachment = Texture | RenderPassColorLayer;

/**
 * The attachments a render pass draws into, named as a whole (as WebGPU's render pass descriptor does): `color` in
 * order (draw buffers), `depth` a depth or depth-stencil texture, or a renderbuffer (its own attachment point). A pass
 * keeps one per target set, created with its textures (`GraphicsResources.renderPassAttachments`), and may change its
 * fields between passes (a ping-pong target); it is a resource, released by `destroy`. How a backend draws into them
 * is its own (WebGL: a framebuffer object each, WebGPU: an attachment set).
 */
export interface RenderPassAttachments {
  readonly kind: 'render-pass-attachments';
  color: readonly RenderPassColorAttachment[];
  depth: Texture | Renderbuffer | undefined;
  destroy(): void;
}

export function isRenderPassAttachments(target: unknown): target is RenderPassAttachments {
  return (
    typeof target === 'object' &&
    target !== null &&
    (target as Partial<RenderPassAttachments>).kind === 'render-pass-attachments'
  );
}

/** Render pass attachments of `resources` (`GraphicsResources.renderPassAttachments`) */
export function createRenderPassAttachments(
  resources: GraphicsResources,
  color: readonly RenderPassColorAttachment[],
  depth?: Texture | Renderbuffer,
): RenderPassAttachments {
  return resources.renderPassAttachments(color, depth);
}

export function isTextureAttachment(t: Texture | Renderbuffer): t is Texture {
  return 'define' in t;
}

/** Calls `f` with each texture and renderbuffer `attachments` names (layers as their texture), allocating nothing */
export function forEachRenderPassAttachment(
  attachments: RenderPassAttachments,
  f: (resource: Texture | Renderbuffer) => void,
) {
  const { color, depth } = attachments;
  for (let i = 0, il = color.length; i < il; ++i) {
    const c = color[i];
    f('layer' in c ? c.texture : c);
  }
  if (depth) f(depth);
}

/**
 * What a backend draws attachments with (its target `T`): how it attaches and detaches them; attaching what a
 * point holds already is skipped by the backends
 */
export interface RenderPassAttachmentsBackend<T> {
  createTarget(): T;
  destroyTarget(target: T): void;
  attachColor(target: T, texture: Texture, index: number, layer: number | undefined): void;
  detachColor(target: T, texture: Texture, index: number): void;
  attachDepth(target: T, depth: Texture | Renderbuffer): void;
  detachDepth(target: T, depth: Texture | Renderbuffer): void;
  setDrawBuffers(target: T, count: number): void;
}

/** Attachments of a backend: `getBackendTarget` (for the backend only) brings its target in line with the fields */
export interface BackendRenderPassAttachments<T> extends RenderPassAttachments {
  getBackendTarget(): T;
}

const MaxColorAttachments = 8;

/** Render pass attachments drawn with a target of `backend`'s of their own (made when first begun) */
export function createBackendRenderPassAttachments<T>(
  backend: RenderPassAttachmentsBackend<T>,
  color: readonly RenderPassColorAttachment[],
  depth?: Texture | Renderbuffer,
): BackendRenderPassAttachments<T> {
  let target: T | undefined;
  /** what the target holds, as last brought in line */
  const attachedColor: (RenderPassColorAttachment | undefined)[] = [];
  let attachedDepth: Texture | Renderbuffer | undefined;
  let drawBuffers = 1;
  let destroyed = false;

  const attachments: BackendRenderPassAttachments<T> = {
    kind: 'render-pass-attachments',
    color,
    depth,
    getBackendTarget: () => {
      if (destroyed) throw new Error('render pass attachments destroyed');
      target ??= backend.createTarget();
      const { color, depth } = attachments;
      if (color.length > MaxColorAttachments) throw new Error(`more than ${MaxColorAttachments} color attachments`);
      for (let i = 0; i < color.length; ++i) {
        const c = color[i];
        if ('layer' in c) backend.attachColor(target, c.texture, i, c.layer);
        else backend.attachColor(target, c, i, undefined);
        attachedColor[i] = c;
      }
      for (let i = color.length; i < attachedColor.length; ++i) {
        const c = attachedColor[i];
        if (c) backend.detachColor(target, 'layer' in c ? c.texture : c, i);
      }
      attachedColor.length = color.length;
      if (attachedDepth && attachedDepth !== depth) backend.detachDepth(target, attachedDepth);
      if (depth) backend.attachDepth(target, depth);
      attachedDepth = depth;
      const count = Math.max(1, color.length);
      if (count !== drawBuffers) {
        backend.setDrawBuffers(target, count);
        drawBuffers = count;
      }
      return target;
    },
    destroy: () => {
      if (destroyed) return;
      destroyed = true;
      if (target !== undefined) backend.destroyTarget(target);
    },
  };
  return attachments;
}

/**
 * Attachments of a target a backend draws otherwise (the drawing buffer's render targets, `createDrawTarget`): never
 * begun as attachments, naming nothing
 */
export function createFixedRenderPassAttachments<T>(target: T): BackendRenderPassAttachments<T> {
  return {
    kind: 'render-pass-attachments',
    color: [],
    depth: undefined,
    getBackendTarget: () => target,
    destroy: () => {},
  };
}
