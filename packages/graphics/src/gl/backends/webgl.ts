/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsBackendDescriptor } from '../backend.js';
import { WebGL2Backend, WebGLBackend } from '../webgl/backend.js';

/** Backend set WebGL: WebGL2, then WebGL */
export function getDefaultGraphicsBackends(): GraphicsBackendDescriptor[] {
  return [WebGL2Backend, WebGLBackend];
}
