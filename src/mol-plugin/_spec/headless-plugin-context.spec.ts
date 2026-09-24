/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import gl from 'gl';
import { HeadlessPluginContext } from '../headless-plugin-context';
import type { ExternalModules } from '../util/headless-screenshot';
import { DefaultPluginSpec } from '../spec';

describe('HeadlessPluginContext', () => {
    it('HeadlessPluginContext', async () => {
        const externalModules: ExternalModules = { gl };
        const spec = DefaultPluginSpec();
        const plugin = new HeadlessPluginContext(externalModules, spec, { height: 100, width: 100 }, {});
        expect(plugin).toBeTruthy();

        await plugin.init();

        plugin.dispose();
    });
});
