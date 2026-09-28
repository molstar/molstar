/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { HeadlessPluginContext } from '../headless-plugin-context';
import type { ExternalModules } from '../util/headless-screenshot';
import { DefaultPluginSpec } from '../spec';

let hasGl = true;
try {
    require.resolve('gl');
} catch (e) {
    if ((e as NodeJS.ErrnoException).code !== 'MODULE_NOT_FOUND') throw e;
    hasGl = false;
    console.warn('Skipping HeadlessPluginContext test: optional package "gl" is not installed.');
}

describe('HeadlessPluginContext', () => {
    (hasGl ? it : it.skip)('HeadlessPluginContext', async () => {
        const gl = require('gl') as typeof import('gl');
        const externalModules: ExternalModules = { gl };
        const spec = DefaultPluginSpec();
        const plugin = new HeadlessPluginContext(externalModules, spec, { height: 100, width: 100 }, {});
        expect(plugin).toBeTruthy();

        await plugin.init();

        plugin.dispose();
    });
});
