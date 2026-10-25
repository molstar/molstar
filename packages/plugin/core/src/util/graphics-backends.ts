/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getDefaultGraphicsBackends as getAllBackends } from '@molstar/graphics/gl/backends/all';
import { getDefaultGraphicsBackends as getWebGLBackends } from '@molstar/graphics/gl/backends/webgl';
import { getDefaultGraphicsBackends as getWebGPUBackends } from '@molstar/graphics/gl/backends/webgpu';
import type { GraphicsBackendDescriptor } from '@molstar/graphics/gl/backend';

/**
 * `auto`: WebGPU where it has the features of WebGL2 (blending into float32 targets, missing in older browsers), else
 * WebGL2, then WebGL1; `webgpu`: WebGPU, else WebGL; `webgl`: WebGL2, else WebGL1; `webgl2` and `webgl1`: only that
 * version
 */
export type GraphicsBackendChoice = 'auto' | 'webgpu' | 'webgl' | 'webgl2' | 'webgl1';

/**
 * The backends of `choice` for `PluginConfig.General.GraphicsBackends`, undefined for any other value (the plugin's default, `auto`);
 * from the backend sets (not the backends), so that builds with one backend alias the others away (see
 * `molstar.variants` in the Viewer's manifest)
 */
export function getGraphicsBackends(choice: GraphicsBackendChoice): readonly GraphicsBackendDescriptor[] | undefined {
  switch (choice) {
    case 'auto':
      return getAllBackends();
    case 'webgl':
      return getWebGLBackends();
    case 'webgpu':
      return [...new Set([...getWebGPUBackends(), ...getWebGLBackends()])];
    case 'webgl2':
    case 'webgl1': {
      const id = choice === 'webgl2' ? 'webgl2' : 'webgl';
      const backends = getWebGLBackends().filter((b) => b.id === id);
      // a build without WebGL has none: its own backends
      return backends.length > 0 ? backends : getAllBackends();
    }
    default:
      return undefined;
  }
}
