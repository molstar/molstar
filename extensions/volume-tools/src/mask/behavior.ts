/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { MaskSelectionColorThemeProvider } from './theme.js';

/** PluginBehavior that marks the mask tool as active in the plugin. */
export const VolumeMaskBehavior = PluginBehavior.create({
  name: 'volume-mask-behavior',
  category: 'misc',
  display: { name: 'Volume Mask Creator' },
  ctor: class extends PluginBehavior.Handler {
    // MaskVolumeFromSource is a BuiltIn transformer, registered at module load.
    private unregisterEntry: (() => void) | undefined;

    register() {
      const entry: PluginRegistryEntry = { volume: { themes: { color: [MaskSelectionColorThemeProvider] } } };
      this.unregisterEntry = this.ctx.register(entry);
    }
    unregister() {
      this.unregisterEntry?.();
      this.unregisterEntry = undefined;
    }
  },
  params: () => ({}),
});
