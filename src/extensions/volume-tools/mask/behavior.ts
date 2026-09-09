/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { PluginBehavior } from '../../../mol-plugin/behavior/behavior';
import { MaskSelectionColorThemeProvider } from './theme';

/** PluginBehavior that marks the mask tool as active in the plugin. */
export const VolumeMaskBehavior = PluginBehavior.create({
    name: 'volume-mask-behavior',
    category: 'misc',
    display: { name: 'Volume Mask Creator' },
    ctor: class extends PluginBehavior.Handler {
        // MaskVolumeFromSource is a BuiltIn transformer, registered at module load.
        register() {
            this.ctx.representation.volume.themes.colorThemeRegistry.add(MaskSelectionColorThemeProvider);
        }
        unregister() {
            this.ctx.representation.volume.themes.colorThemeRegistry.remove(MaskSelectionColorThemeProvider);
        }
    },
    params: () => ({}),
});
