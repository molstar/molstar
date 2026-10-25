/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { AssetManager } from '@molstar/core/util/assets';
import { createContext, getGLContext } from '@molstar/graphics/gl/webgl/context';
import { WebGL2Backend, WebGLBackend } from '@molstar/graphics/gl/webgl/backend';
import { Canvas3DContext } from './canvas3d.js';

/**
 * A `Canvas3DContext` with a WebGL context created synchronously on `canvas` (WebGL2 unless unavailable); in its own
 * module, so that bundles importing `Canvas3D` without it contain no WebGL code. For other
 * backends, or a choice among several, see `Canvas3DContext.create`.
 */
export function createWebGLCanvas3DContext(
  canvas: HTMLCanvasElement,
  assetManager: AssetManager,
  attribs: Partial<Canvas3DContext.Attribs> = {},
  props: Partial<Canvas3DContext.Props> = {},
): Canvas3DContext {
  const a = { ...Canvas3DContext.DefaultAttribs, ...attribs };
  const p = { ...Canvas3DContext.DefaultProps, ...props };
  const { powerPreference, failIfMajorPerformanceCaveat, antialias, preserveDrawingBuffer } = a;
  const gl = getGLContext(canvas, {
    powerPreference,
    failIfMajorPerformanceCaveat,
    antialias,
    preserveDrawingBuffer,
    alpha: true, // the renderer requires an alpha channel
    depth: true, // the renderer requires a depth buffer
    stencil: true, // the renderer requires a stencil buffer
    premultipliedAlpha: true, // the renderer outputs PMA
  });
  if (gl === null) throw new Error('Could not create a WebGL rendering context');
  const graphics = createContext(gl, { pixelScale: Canvas3DContext.getPixelScale(p) });
  // to create a context again after a loss (restore, fallback to WebGL1)
  return Canvas3DContext.fromGraphics(canvas, graphics, assetManager, a, p, [WebGL2Backend, WebGLBackend]);
}
