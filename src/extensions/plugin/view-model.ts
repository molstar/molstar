/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultPluginSpec, PluginSpec } from '../../mol-plugin/spec';
import { PluginContext } from '../../mol-plugin/context';
import { SingleAsyncQueue } from '../../mol-util/single-async-queue';

export class PluginViewModel {
    private mountQueue = new SingleAsyncQueue();
    readonly plugin: PluginContext;

    get initialized() {
        return this.plugin.initialized;
    }

    private async init() {
        await this.plugin.init();
    }

    mount(root: HTMLElement) {
        this.mountQueue.enqueue(() => this.plugin.mountAsync(root));
    }

    unmount() {
        this.mountQueue.enqueue(() => this.plugin.unmount());
    }

    constructor(options?: { spec?: PluginSpec }) {
        const spec = options?.spec ?? DefaultPluginSpec();
        this.plugin = new PluginContext(spec);
        this.init();
    }
}