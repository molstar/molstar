/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsBackendDescriptor } from '../backend.js';
import { WebGL2Backend, WebGLBackend } from '../webgl/backend.js';
import { AutoWebGPUBackend } from '../webgpu/backend.js';

/**
 * Backend set all: WebGPU where it has the features of WebGL2 (`AutoWebGPUBackend`), else WebGL2, then WebGL
 */
export function getDefaultGraphicsBackends(): GraphicsBackendDescriptor[] {
  return [AutoWebGPUBackend, WebGL2Backend, WebGLBackend];
}
