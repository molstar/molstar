/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRequire } from 'node:module';
import { pathToFileURL } from 'node:url';
import type { ExternalModules } from '../screenshot.js';

export type WebGPUModule = NonNullable<ExternalModules['webgpu']>;

/**
 * Dawn's Node.js bindings, the optional `webgpu` peer: installed with the workspace, or by
 * `pnpm native:install -- --webgpu` and found through `NODE_PATH` (`pnpm test:native`, `pnpm native:run`), which
 * `import` ignores, so it is resolved with `require.resolve`
 */
async function importWebGPU(): Promise<WebGPUModule> {
  try {
    return (await import('webgpu' as string)) as WebGPUModule;
  } catch {
    const path = createRequire(import.meta.url).resolve('webgpu');
    return (await import(pathToFileURL(path).href)) as WebGPUModule;
  }
}

/**
 * The `webgpu` module with a usable adapter, or why there is none: tests skip with the reason
 * (`describe.skipIf(!WebGPU.module)`), so that they show as skipped instead of passing
 */
export const WebGPU: { module?: WebGPUModule; skipReason?: string } = await (async () => {
  let module: WebGPUModule;
  try {
    module = await importWebGPU();
  } catch {
    return { skipReason: 'the optional package "webgpu" is not installed (pnpm native:install -- --webgpu)' };
  }
  const adapter = await module.create([]).requestAdapter();
  if (!adapter) return { skipReason: 'no WebGPU adapter' };
  return { module };
})();

if (WebGPU.skipReason) console.warn(`Skipping WebGPU tests: ${WebGPU.skipReason}.`);
