/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import type { ExternalModules } from '@molstar/plugin-headless/screenshot';
import { DefaultPluginSpec } from '@molstar/plugin/spec';

let hasGl = true;
let glModule: typeof import('gl') | undefined;
try {
    glModule = require('gl') as typeof import('gl');
    const context = glModule(1, 1);
    if (!context) throw new Error('Optional package "gl" loaded but could not create a WebGL context.');
    context.getExtension('STACKGL_destroy_context')?.destroy();
} catch (e) {
    const error = e as NodeJS.ErrnoException;
    const message = error.message ?? String(e);
    const packageMissing = error.code === 'MODULE_NOT_FOUND' && /Cannot find module ['"]gl['"]/.test(message);
    const bindingMissing = /Could not locate the bindings file/i.test(message)
        || (error.code === 'ERR_DLOPEN_FAILED' && /NODE_MODULE_VERSION|compiled against|different Node\.js version/i.test(message))
        || (error.code === 'MODULE_NOT_FOUND' && /\.node['"]/.test(message) && (error.requireStack ?? []).some(path => /[/\\]gl[/\\]/.test(path)));
    if (!packageMissing && !bindingMissing) throw e;
    hasGl = false;
    console.warn(`Skipping HeadlessPluginContext test: optional native package "gl" is unavailable (${message}).`);
}

describe('HeadlessPluginContext', () => {
    (hasGl ? it : it.skip)('HeadlessPluginContext', async () => {
        const externalModules: ExternalModules = { gl: glModule! };
        const spec = DefaultPluginSpec();
        const plugin = new HeadlessPluginContext(externalModules, spec, { height: 100, width: 100 }, {});
        expect(plugin).toBeTruthy();

        await plugin.init();

        plugin.dispose();
    });
});
