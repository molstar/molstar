/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * The backends used when none are given (e.g. by the plugin without `PluginConfig.General.GraphicsBackends`): the `all`
 * set (WebGPU where capable, else WebGL2, then WebGL). Builds with one backend alias this module to that backend's set,
 * e.g. with esbuild `alias: { '@molstar/graphics/gl/backends/default': '@molstar/graphics/gl/backends/webgpu' }`, so
 * that the WebGL backend is not bundled.
 */
export { getDefaultGraphicsBackends } from './all.js';
