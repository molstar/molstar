/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Framebuffer } from './framebuffer.js';

/** depth precision is a preference, backends fall back to the best supported depth format */
export type RenderbufferFormat =
  | 'depth16'
  | 'stencil8'
  | 'rgba4'
  | 'depth-stencil'
  | 'depth24'
  | 'depth32f'
  | 'depth24-stencil8'
  | 'depth32f-stencil8';
export type RenderbufferAttachment = 'depth' | 'stencil' | 'depth-stencil' | 'color0';

export interface Renderbuffer {
  readonly id: number;

  getByteCount: () => number;

  bind: () => void;
  attachFramebuffer: (framebuffer: Framebuffer) => void;
  detachFramebuffer: (framebuffer: Framebuffer) => void;
  setSize: (width: number, height: number) => void;
  /** Frees the GPU storage, keeping the size: the next use allocates it again, see `Texture.releaseStorage` */
  releaseStorage: () => void;
  destroy: () => void;
}
