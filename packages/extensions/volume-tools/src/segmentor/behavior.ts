/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { VolumeSegmentorManager } from './manager.js';
import { BodyLabelColorThemeProvider } from './theme.js';

/**
 * Registers the body-label color theme and creates the `VolumeSegmentorManager` for the plugin
 * (available via `VolumeSegmentorManager.get(plugin)`). `BodyMaskFromLabels` is a BuiltIn
 * transformer, registered at module load. For a full interactive UI, see `src/examples/volume-tools/`.
 */
export const VolumeSegmentorBehavior = PluginBehavior.create({
  name: 'volume-segmentor',
  category: 'misc',
  display: {
    name: 'Volume Segmentor',
    description: 'Interactive segmentation of a volume into soft-edged body masks.',
  },
  ctor: class extends PluginBehavior.Handler {
    private manager: VolumeSegmentorManager | undefined;
    private unregisterEntry: (() => void) | undefined;

    register() {
      this.manager = new VolumeSegmentorManager(this.ctx);
      VolumeSegmentorManager.register(this.ctx, this.manager);
      const entry: PluginRegistryEntry = { volume: { themes: { color: [BodyLabelColorThemeProvider] } } };
      this.unregisterEntry = this.ctx.register(entry);
    }

    unregister() {
      this.unregisterEntry?.();
      this.unregisterEntry = undefined;
      VolumeSegmentorManager.unregister(this.ctx);
      this.manager?.dispose();
      this.manager = undefined;
    }
  },
  params: () => ({}),
});
