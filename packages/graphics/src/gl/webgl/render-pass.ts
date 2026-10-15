/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { isDebugMode } from '@molstar/core/util/debug';
import type { WebGLState } from './state.js';

/**
 * Ends the open render pass, for operations that can't be part of one: copies, readback and
 * attachment changes. On WebGL they also bind another framebuffer, so draws need a new pass.
 */
export function endRenderPass(state: WebGLState) {
  if (state.inRenderPass) state.onRenderPassEnd?.();
  state.inRenderPass = false;
}

/** Throws in debug mode when drawing or clearing outside of a render pass */
export function checkRenderPass(state: WebGLState, operation: string) {
  if (isDebugMode && !state.inRenderPass) {
    throw new Error(`${operation} outside of a render pass, call 'beginRenderPass' first`);
  }
}
