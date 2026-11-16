/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import type { TextureAttachment } from '../texture.js';

/** A texture or renderbuffer by its id (the resources may wrap the objects they create, ids stay) */
export interface AttachmentRef {
  readonly kind: 'texture' | 'renderbuffer';
  readonly id: number;
}

const getNextAttachmentSetId = idFactory();

/** A texture or renderbuffer in an attachment set */
export interface WebGPUAttachment {
  readonly format: GPUTextureFormat;
  readonly hasStencil: boolean;
  /** the slice of a 3d texture rendered into */
  readonly depthSlice?: number;
  getAttachmentView: () => GPUTextureView | undefined;
  getGPUTexture: () => GPUTexture | undefined;
  getWidth: () => number;
  getHeight: () => number;
  /** what it is of, for `WebGPUAttachmentSet.getAttachments` (layer attachments: their texture's) */
  readonly ref: AttachmentRef;
}

/**
 * The attachments of a render pass target, kept as a set the render pass descriptor is made from (WebGPU has no
 * framebuffer objects: passes name their attachments)
 */
export interface WebGPUAttachmentSet {
  readonly id: number;
  /** Number of color attachments written by draws; default 1 */
  setDrawBuffers: (count: number) => void;
  /** The textures and renderbuffers attached (layers of 3d textures as their texture) */
  getAttachments: () => readonly AttachmentRef[];
  destroy: () => void;
  attach: (attachment: TextureAttachment, target: WebGPUAttachment) => void;
  detach: (attachment: TextureAttachment) => void;
  /** color attachments `0` to `drawBuffers - 1`, empty for depth-only targets */
  getColorAttachments: () => (WebGPUAttachment | undefined)[];
  getColorAttachment: (index: number) => WebGPUAttachment | undefined;
  getDepthAttachment: () => WebGPUAttachment | undefined;
  /** changes with every attach, detach and draw buffer count: a pass target made of these attachments stays valid */
  readonly version: number;
}

const AttachmentSets = new WeakSet<object>();

export function isWebGPUAttachmentSet(set: object): set is WebGPUAttachmentSet {
  return AttachmentSets.has(set);
}

/** Narrows to `WebGPUAttachmentSet`, throws for objects of another backend */
export function asWebGPUAttachmentSet(set: object): WebGPUAttachmentSet {
  if (!isWebGPUAttachmentSet(set)) throw new Error('expected a WebGPU attachment set');
  return set;
}

function getColorIndex(attachment: TextureAttachment): number | undefined {
  if (typeof attachment === 'number') return attachment;
  // 'color0' to 'color7', without a regular expression (attached every pass: no match objects)
  if (attachment.length === 6 && attachment.startsWith('color')) return attachment.charCodeAt(5) - 48;
  return undefined;
}

export function createWebGPUAttachmentSet(): WebGPUAttachmentSet {
  const colors: (WebGPUAttachment | undefined)[] = [];
  let depth: WebGPUAttachment | undefined;
  let drawBuffers = 1;
  let version = 0;

  const set: WebGPUAttachmentSet = {
    id: getNextAttachmentSetId(),
    setDrawBuffers: (count: number) => {
      if (count !== drawBuffers) version += 1;
      drawBuffers = count;
    },
    attach: (attachment, texture) => {
      const index = getColorIndex(attachment);
      const current = index !== undefined ? colors[index] : depth;
      // the same texture (slice) again: no change
      if (
        current &&
        current.ref.kind === texture.ref.kind &&
        current.ref.id === texture.ref.id &&
        current.depthSlice === texture.depthSlice
      ) {
        if (current !== texture) {
          if (index !== undefined) colors[index] = texture;
          else depth = texture;
        }
        return;
      }
      if (index !== undefined) colors[index] = texture;
      else depth = texture;
      version += 1;
    },
    detach: (attachment) => {
      const index = getColorIndex(attachment);
      if ((index !== undefined ? colors[index] : depth) === undefined) return;
      if (index !== undefined) colors[index] = undefined;
      else depth = undefined;
      version += 1;
    },
    get version() {
      return version;
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
  AttachmentSets.add(set);
  return set;
}
