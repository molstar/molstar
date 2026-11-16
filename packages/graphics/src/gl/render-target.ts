/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import type { GraphicsResources, RenderTargetDepthStencil, RenderTargetFormat, RenderTargetType } from './context.js';
import type { Renderbuffer } from './renderbuffer.js';
import type { Texture, TextureFilter } from './texture.js';
import { createRenderPassAttachments, type RenderPassAttachments } from './render-pass-attachments.js';

export interface RenderTarget {
  readonly id: number;
  readonly texture: Texture;
  readonly depthRenderbuffer: Renderbuffer | null;
  /** the texture and the depth renderbuffer, as a pass draws into them (`beginRenderPass` with the target) */
  readonly attachments: RenderPassAttachments;

  getByteCount: () => number;

  getWidth: () => number;
  getHeight: () => number;
  setSize: (width: number, height: number) => void;
  destroy: () => void;
}

const getNextRenderTargetId = idFactory();

/** A color texture with an optional depth (and stencil) renderbuffer, from the resources of any backend */
export function createRenderTargetFromResources(
  resources: GraphicsResources,
  _width: number,
  _height: number,
  depthStencil: RenderTargetDepthStencil = 'depth',
  type: RenderTargetType = 'uint8',
  filter: TextureFilter = 'nearest',
  format: RenderTargetFormat = 'rgba',
): RenderTarget {
  const targetTexture =
    type === 'fp16'
      ? resources.texture('image-float16', format, 'fp16', filter)
      : type === 'float32'
        ? resources.texture('image-float32', format, 'float', filter)
        : resources.texture('image-uint8', format, 'ubyte', filter);
  // make a depth renderbuffer of the same size as the targetTexture
  const depthRenderbuffer =
    depthStencil === 'none'
      ? null
      : depthStencil === 'depth-stencil'
        ? resources.renderbuffer('depth32f-stencil8', 'depth-stencil', _width, _height)
        : resources.renderbuffer('depth32f', 'depth', _width, _height);

  targetTexture.define(_width, _height);
  const attachments = createRenderPassAttachments(resources, [targetTexture], depthRenderbuffer ?? undefined);

  let destroyed = false;

  return {
    id: getNextRenderTargetId(),
    texture: targetTexture,
    depthRenderbuffer,
    attachments,

    getByteCount: () => targetTexture.getByteCount() + (depthRenderbuffer ? depthRenderbuffer.getByteCount() : 0),

    getWidth: () => _width,
    getHeight: () => _height,
    setSize: (width: number, height: number) => {
      if (_width === width && _height === height) {
        return;
      }

      _width = width;
      _height = height;
      targetTexture.define(_width, _height);
      if (depthRenderbuffer) depthRenderbuffer.setSize(_width, _height);
    },
    destroy: () => {
      if (destroyed) return;
      targetTexture.destroy();
      attachments.destroy();
      if (depthRenderbuffer) depthRenderbuffer.destroy();
      destroyed = true;
    },
  };
}
