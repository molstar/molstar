/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { PluginBehavior } from '../../../mol-plugin/behavior/behavior';
import { VolumeSegmentorManager } from './manager';
import { BodyLabelColorThemeProvider } from './theme';

/**
 * Registers the body-label color theme and creates the `VolumeSegmentorManager` for the plugin
 * (available via `VolumeSegmentorManager.get(plugin)`). `BodyMaskFromLabels` is a BuiltIn
 * transformer, registered at module load. For a full interactive UI, see `src/examples/volume-tools/`.
 */
export const VolumeSegmentorBehavior = PluginBehavior.create({
    name: 'volume-segmentor',
    category: 'misc',
    display: { name: 'Volume Segmentor', description: 'Interactive segmentation of a volume into soft-edged body masks.' },
    ctor: class extends PluginBehavior.Handler {
        private manager: VolumeSegmentorManager | undefined;

        register() {
            this.manager = new VolumeSegmentorManager(this.ctx);
            VolumeSegmentorManager.register(this.ctx, this.manager);
            this.ctx.representation.volume.themes.colorThemeRegistry.add(BodyLabelColorThemeProvider);
        }

        unregister() {
            this.ctx.representation.volume.themes.colorThemeRegistry.remove(BodyLabelColorThemeProvider);
            VolumeSegmentorManager.unregister(this.ctx);
            this.manager?.dispose();
            this.manager = undefined;
        }
    },
    params: () => ({}),
});
