/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { setCanvasModule } from '../../mol-geo/geometry/text/font-atlas';
import { ExternalModules } from './headless-screenshot';

/** Optional Node dependencies are loaded only by headless applications. */
export async function loadHeadlessModules(backend: 'webgpu' | 'webgl' = 'webgpu'): Promise<ExternalModules> {
    const modules: ExternalModules = {};
    if (backend === 'webgpu') {
        // Keep the Node binding's ambient @webgpu/types out of the browser DOM build.
        const moduleName: string = 'webgpu';
        const { create, globals } = await import(moduleName) as { create(options: string[]): GPU, globals: Record<string, unknown> };
        Object.assign(globalThis, globals);
        modules.webgpu = create([]);
    } else {
        modules.gl = (await import('gl')).default;
    }
    modules.pngjs = await import('pngjs');
    modules['jpeg-js'] = await import('jpeg-js');
    const canvasModule: string = '@napi-rs/canvas';
    setCanvasModule(await import(canvasModule));
    return modules;
}
