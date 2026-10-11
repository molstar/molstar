/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { TextureFilter } from '../texture.js';
import { createRenderTargetFromResources, type RenderTarget } from '../render-target.js';
import type { WebGLResources } from './resources.js';
import { type GLRenderingContext, isWebGL2 } from './compat.js';

export function createRenderTarget(
  gl: GLRenderingContext,
  resources: WebGLResources,
  width: number,
  height: number,
  depthStencil: 'none' | 'depth' | 'depth-stencil' = 'depth',
  type: 'uint8' | 'float32' | 'fp16' = 'uint8',
  filter: TextureFilter = 'nearest',
  format: 'rgba' | 'r' = 'rgba',
): RenderTarget {
  // WebGL1 can't render to `r` (luminance), `rgba` keeps the value in `.r` as well
  if (format === 'r' && !isWebGL2(gl)) format = 'rgba';
  return createRenderTargetFromResources(resources, width, height, depthStencil, type, filter, format);
}
