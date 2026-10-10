/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginUISpec } from '@molstar/plugin-ui/spec';
import { PluginUIContext } from '@molstar/plugin-ui/context';

export class PluginUIViewModel {
  readonly plugin: PluginUIContext;

  get initialized() {
    return this.plugin.initialized;
  }

  private async init() {
    await this.plugin.init();
  }

  /** The spec is required: pass `DefaultPluginUISpec()` for the default plugin. */
  constructor(options: { spec: PluginUISpec }) {
    this.plugin = new PluginUIContext(options.spec);
    this.init();
  }
}
