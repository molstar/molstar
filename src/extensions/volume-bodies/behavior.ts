/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { VolumeBodiesManager } from './manager';
import { BodyLabelColorThemeProvider } from './theme';

/**
 * Registers the body-label color theme and creates the `VolumeBodiesManager` for the plugin
 * (available via `VolumeBodiesManager.get(plugin)`). `BodyMaskFromLabels` is a BuiltIn
 * transformer, registered at module load. For a full interactive UI, see `src/examples/volume-bodies/`.
 */
export const VolumeBodiesBehavior = PluginBehavior.create({
    name: 'volume-bodies',
    category: 'misc',
    display: { name: 'Volume Bodies', description: 'Interactive segmentation of a volume into soft-edged body masks.' },
    ctor: class extends PluginBehavior.Handler {
        private manager: VolumeBodiesManager | undefined;

        register() {
            this.manager = new VolumeBodiesManager(this.ctx);
            VolumeBodiesManager.register(this.ctx, this.manager);
            this.ctx.representation.volume.themes.colorThemeRegistry.add(BodyLabelColorThemeProvider);
        }

        unregister() {
            this.ctx.representation.volume.themes.colorThemeRegistry.remove(BodyLabelColorThemeProvider);
            VolumeBodiesManager.unregister(this.ctx);
            this.manager?.dispose();
            this.manager = undefined;
        }
    },
    params: () => ({}),
});
