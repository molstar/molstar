/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginSpec } from '@molstar/plugin/spec';
import { PluginContext } from '@molstar/plugin/context';
import { SingleAsyncQueue } from '@molstar/core/util/single-async-queue';

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

  /** The spec is required: pass `DefaultPluginSpec()` for the default plugin. */
  constructor(options: { spec: PluginSpec }) {
    this.plugin = new PluginContext(options.spec);
    this.init();
  }
}
