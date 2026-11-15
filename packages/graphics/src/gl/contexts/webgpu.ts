/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// WebGPU graphics contexts, made and recognized outside the backend: the backend folder (`gl/webgpu/`) is internal to
// the package (its exports source-only), this module the public way in

export {
  createWebGPUContext,
  isWebGPUContext,
  requestWebGPUDevice,
  type WebGPUContext,
  type WebGPUContextOptions,
  type WebGPUFeatures,
  type WebGPULimits,
} from '../webgpu/context.js';
export { AutoWebGPUBackend, createWebGPUBackend, WebGPUBackend } from '../webgpu/backend.js';
export { WebGPUUnsupportedError } from '../webgpu/errors.js';
