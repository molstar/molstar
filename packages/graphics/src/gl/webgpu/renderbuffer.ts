/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import type { Renderbuffer, RenderbufferAttachment, RenderbufferFormat } from '../renderbuffer.js';
import { TextureUsage } from './constants.js';
import { asWebGPUFramebuffer, type WebGPUAttachment } from './framebuffer.js';
import { countResource, createByteAccount, type WebGPUResourceEnv } from './texture.js';

const getNextRenderbufferId = idFactory();

/** Depth precision is a preference: falls back to the formats of WebGPU core */
function getRenderbufferFormat(device: GPUDevice, format: RenderbufferFormat): GPUTextureFormat {
  switch (format) {
    case 'depth16':
      return 'depth16unorm';
    case 'depth24':
      return 'depth24plus';
    case 'depth32f':
      return 'depth32float';
    case 'stencil8':
      return 'stencil8';
    case 'rgba4':
      return 'rgba8unorm';
    case 'depth-stencil':
    case 'depth24-stencil8':
      return 'depth24plus-stencil8';
    case 'depth32f-stencil8':
      return device.features.has('depth32float-stencil8') ? 'depth32float-stencil8' : 'depth24plus-stencil8';
  }
}

function getBytesPerTexel(format: GPUTextureFormat) {
  switch (format) {
    case 'depth16unorm':
      return 2;
    case 'stencil8':
      return 1;
    case 'depth32float-stencil8':
      return 8;
    default:
      return 4;
  }
}

export interface WebGPURenderbuffer extends Renderbuffer, WebGPUAttachment {}

/** Attachment-only texture (not sampled) */
export function createWebGPURenderbuffer(
  env: WebGPUResourceEnv,
  format: RenderbufferFormat,
  attachment: RenderbufferAttachment,
  _width: number,
  _height: number,
): WebGPURenderbuffer {
  const { release } = env;
  const gpuFormat = getRenderbufferFormat(env.device, format);
  const hasStencil = gpuFormat.includes('stencil');

  let width = _width;
  let height = _height;
  let texture: GPUTexture | undefined;
  let view: GPUTextureView | undefined;
  const account = createByteAccount(env.byteCounts, 'renderbuffer');
  const uncount = countResource(env.stats, 'renderbuffer');

  let destroyed = false;

  /** Allocates the storage when first attached to a render pass (lazy: unused attachments take no memory) */
  function ensureStorage() {
    if (texture || destroyed) return;
    account(Math.max(1, width) * Math.max(1, height) * getBytesPerTexel(gpuFormat));
    texture = env.device.createTexture({
      size: { width: Math.max(1, width), height: Math.max(1, height) },
      format: gpuFormat,
      usage: TextureUsage.RENDER_ATTACHMENT,
      label: `renderbuffer ${format}`,
    });
    view = texture.createView();
  }

  /** Releases the storage of the previous size; the next attachment allocates it at the new one */
  function releaseStorage() {
    if (!texture) return;
    release(texture);
    texture = undefined;
    view = undefined;
    account(0);
  }

  const id = getNextRenderbufferId();
  const renderbuffer: WebGPURenderbuffer = {
    id,
    ref: { kind: 'renderbuffer', id },
    format: gpuFormat,
    hasStencil,
    getAttachmentView: () => {
      ensureStorage();
      return view;
    },
    getGPUTexture: () => {
      ensureStorage();
      return texture;
    },
    getWidth: () => width,
    getHeight: () => height,

    getByteCount: () => (texture ? width * height * getBytesPerTexel(gpuFormat) : 0),
    bind: () => {},
    attachFramebuffer: (framebuffer) => {
      asWebGPUFramebuffer(framebuffer).attach(attachment, renderbuffer);
    },
    detachFramebuffer: (framebuffer) => {
      asWebGPUFramebuffer(framebuffer).detach(attachment);
    },
    releaseStorage,
    setSize: (w: number, h: number) => {
      if (w === width && h === height) return;
      width = w;
      height = h;
      releaseStorage();
    },
    destroy: () => {
      if (destroyed) return;
      destroyed = true;
      releaseStorage();
      uncount();
    },
  };
  return renderbuffer;
}
