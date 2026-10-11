/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import type { Renderbuffer, RenderbufferAttachment, RenderbufferFormat } from '../renderbuffer.js';
import { TextureUsage } from './constants.js';
import { asWebGPUFramebuffer, type WebGPUAttachment } from './framebuffer.js';
import type { WebGPUResourceEnv } from './texture.js';

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

  function init() {
    if (texture) release(texture);
    texture = env.device.createTexture({
      size: { width: Math.max(1, width), height: Math.max(1, height) },
      format: gpuFormat,
      usage: TextureUsage.RENDER_ATTACHMENT,
      label: `renderbuffer ${format}`,
    });
    view = texture.createView();
  }
  init();

  const renderbuffer: WebGPURenderbuffer = {
    id: getNextRenderbufferId(),
    format: gpuFormat,
    hasStencil,
    getAttachmentView: () => view,
    getGPUTexture: () => texture,
    getWidth: () => width,
    getHeight: () => height,

    getByteCount: () => width * height * getBytesPerTexel(gpuFormat),
    bind: () => {},
    attachFramebuffer: (framebuffer) => {
      asWebGPUFramebuffer(framebuffer).attach(attachment, renderbuffer);
    },
    detachFramebuffer: (framebuffer) => {
      asWebGPUFramebuffer(framebuffer).detach(attachment);
    },
    setSize: (w: number, h: number) => {
      if (w === width && h === height) return;
      width = w;
      height = h;
      init();
    },
    reset: () => init(),
    destroy: () => {
      if (texture) release(texture);
      texture = undefined;
      view = undefined;
      env.restorables.delete(renderbuffer);
    },
  };
  env.restorables.add(renderbuffer);
  return renderbuffer;
}
