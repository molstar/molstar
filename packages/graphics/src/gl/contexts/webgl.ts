/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// WebGL graphics contexts, made and recognized outside the backend: the backend folder (`gl/webgl/`) is internal to the
// package (its exports source-only), this module the public way in

export {
  createContext,
  getGLContext,
  isWebGLContext,
  type WebGLContext,
  type WebGLContextProps,
} from '../webgl/context.js';
export { createWebGLBackend, WebGL2Backend, WebGLBackend } from '../webgl/backend.js';
