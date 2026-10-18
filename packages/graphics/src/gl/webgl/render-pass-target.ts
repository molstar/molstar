/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '../context.js';

/**
 * Rectangles have a top-left origin in the neutral API, GL's a bottom-left one: the backend flips them with the height
 * of the render pass target, which each context registers here
 */
const TargetHeights = new WeakMap<GraphicsContext, () => number>();

export function registerRenderPassTargetHeight(ctx: GraphicsContext, getHeight: () => number) {
  TargetHeights.set(ctx, getHeight);
}

/** The height of the target of the open (or last) render pass of `ctx` */
export function getRenderPassTargetHeight(ctx: GraphicsContext) {
  const getHeight = TargetHeights.get(ctx);
  if (!getHeight) throw new Error('not a WebGL context');
  return getHeight();
}

/**
 * Viewport in texel rows from the bottom (GL's origin), for the data textures of the backend's compute; call it in
 * the render pass (flipped with the height of its target)
 */
export function setTexelViewport(ctx: GraphicsContext, x: number, y: number, width: number, height: number) {
  ctx.setViewport(x, getRenderPassTargetHeight(ctx) - y - height, width, height);
}

/** Scissor in texel rows from the bottom, see `setTexelViewport` */
export function setTexelScissor(ctx: GraphicsContext, x: number, y: number, width: number, height: number) {
  ctx.setScissor(x, getRenderPassTargetHeight(ctx) - y - height, width, height);
}
