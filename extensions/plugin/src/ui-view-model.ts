/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultPluginUISpec } from '@molstar/plugin-ui/default-spec';
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

  constructor(options?: { spec?: PluginUISpec }) {
    const spec = options?.spec ?? DefaultPluginUISpec();
    this.plugin = new PluginUIContext(spec);
    this.init();
  }
}
