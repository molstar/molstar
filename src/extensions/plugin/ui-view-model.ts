/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultPluginUISpec, PluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginUIContext } from '../../mol-plugin-ui/context';

export class PluginUIViewModel {
    readonly plugin: PluginUIContext;

    get initialized() {
        return this.plugin.initialized;
    }

    private async init() {
        await this.plugin.init();
    }

    constructor(options?: { spec?: PluginUISpec }) {
        const spec = options?.spec ?? DefaultPluginUISpec();
        this.plugin = new PluginUIContext(spec);
        this.init();
    }
}