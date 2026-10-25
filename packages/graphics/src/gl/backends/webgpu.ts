/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsBackendDescriptor } from '../backend.js';
import { WebGPUBackend } from '../webgpu/backend.js';

/** Backend set WebGPU: WebGPU only, for builds without WebGL (see `./default`) */
export function getDefaultGraphicsBackends(): GraphicsBackendDescriptor[] {
  return [WebGPUBackend];
}
