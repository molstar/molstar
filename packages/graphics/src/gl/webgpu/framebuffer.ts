/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import type { AttachmentRef, Framebuffer } from '../framebuffer.js';
import type { TextureAttachment } from '../texture.js';

const getNextFramebufferId = idFactory();

/** A texture or renderbuffer attached to a framebuffer */
export interface WebGPUAttachment {
  readonly format: GPUTextureFormat;
  readonly hasStencil: boolean;
  /** the slice of a 3d texture rendered into */
  readonly depthSlice?: number;
  getAttachmentView: () => GPUTextureView | undefined;
  getGPUTexture: () => GPUTexture | undefined;
  getWidth: () => number;
  getHeight: () => number;
  /** what it is of, for `Framebuffer.getAttachments` (layer attachments: their texture's) */
  readonly ref: AttachmentRef;
}

/** Attachments of a render pass target; WebGPU has no framebuffer objects, passes name their attachments */
export interface WebGPUFramebuffer extends Framebuffer {
  attach: (attachment: TextureAttachment, target: WebGPUAttachment) => void;
  detach: (attachment: TextureAttachment) => void;
  /** color attachments `0` to `drawBuffers - 1`, empty for depth-only targets */
  getColorAttachments: () => (WebGPUAttachment | undefined)[];
  getColorAttachment: (index: number) => WebGPUAttachment | undefined;
  getDepthAttachment: () => WebGPUAttachment | undefined;
}

const Framebuffers = new WeakSet<Framebuffer>();

export function isWebGPUFramebuffer(framebuffer: Framebuffer): framebuffer is WebGPUFramebuffer {
  return Framebuffers.has(framebuffer);
}

/** Narrows to `WebGPUFramebuffer`, throws for framebuffers of another backend */
export function asWebGPUFramebuffer(framebuffer: Framebuffer): WebGPUFramebuffer {
  if (!isWebGPUFramebuffer(framebuffer)) throw new Error('expected a WebGPU framebuffer');
  return framebuffer;
}

function getColorIndex(attachment: TextureAttachment): number | undefined {
  if (typeof attachment === 'number') return attachment;
  const m = attachment.match(/^color(\d)$/);
  return m ? parseInt(m[1]) : undefined;
}

export function createWebGPUFramebuffer(): WebGPUFramebuffer {
  const colors: (WebGPUAttachment | undefined)[] = [];
  let depth: WebGPUAttachment | undefined;
  let drawBuffers = 1;

  const framebuffer: WebGPUFramebuffer = {
    id: getNextFramebufferId(),
    bind: () => {},
    setDrawBuffers: (count: number) => {
      drawBuffers = count;
    },
    attach: (attachment, texture) => {
      const index = getColorIndex(attachment);
      if (index !== undefined) colors[index] = texture;
      else depth = texture;
    },
    detach: (attachment) => {
      const index = getColorIndex(attachment);
      if (index !== undefined) colors[index] = undefined;
      else depth = undefined;
    },
    getColorAttachments: () => {
      const out: (WebGPUAttachment | undefined)[] = [];
      for (let i = 0; i < drawBuffers; ++i) out.push(colors[i]);
      // a depth-only target has no color attachment
      return out.some((t) => t) ? out : [];
    },
    getColorAttachment: (index) => colors[index],
    getAttachments: () => {
      const out: AttachmentRef[] = [];
      for (const a of colors) if (a) out.push(a.ref);
      if (depth) out.push(depth.ref);
      return out;
    },
    getDepthAttachment: () => depth,
    destroy: () => {
      colors.length = 0;
      depth = undefined;
    },
  };
  Framebuffers.add(framebuffer);
  return framebuffer;
}
