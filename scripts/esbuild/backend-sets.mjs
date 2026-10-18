/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// Builds with a subset of the graphics backends: the backend set modules (`@molstar/graphics/gl/backends/*`, each
// exporting `getDefaultGraphicsBackends`) other than the kept one are aliased to it, so that code choosing among the
// backends (the plugin default, the Viewer's `graphics-backend` option) bundles only the kept backend.

const BackendSetModules = ['default', 'webgl', 'webgpu', 'all'];

/**
 * esbuild `alias` keeping only the backend set `set` (`webgl` or `webgpu`). esbuild resolves alias targets from the
 * working directory: `target` maps the kept module to a path that resolves there (by default its package path)
 */
export function getBackendSetAlias(set, target = (module) => `@molstar/graphics/gl/backends/${module}`) {
  if (!['webgl', 'webgpu'].includes(set)) throw new Error(`unknown backend set '${set}', expected webgl or webgpu`);
  const alias = {};
  for (const module of BackendSetModules) {
    if (module !== set) alias[`@molstar/graphics/gl/backends/${module}`] = target(set);
  }
  return alias;
}
