/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ConfalPyramidsPreset } from './confal-pyramids/behavior';
import { ConfalPyramidsColorThemeProvider } from './confal-pyramids/color';
import { ConfalPyramidsProvider } from './confal-pyramids/property';
import { ConfalPyramidsRepresentationProvider } from './confal-pyramids/representation';
import { NtCTubePreset } from './ntc-tube/behavior';
import { NtCTubeColorThemeProvider } from './ntc-tube/color';
import { NtCTubeProvider } from './ntc-tube/property';
import { NtCTubeRepresentationProvider } from './ntc-tube/representation';


export const DnatcoNtCs = PluginBehavior.create<{ autoAttach: boolean, showToolTip: boolean }>({
    name: 'dnatco-ntcs',
    category: 'custom-props',
    display: {
        name: 'DNATCO NtC Annotations',
        description: 'DNATCO NtC Annotations',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showToolTip: boolean }> {
        register(): void {
            this.ctx.customModelProperties.register(ConfalPyramidsProvider, this.params.autoAttach);
            this.ctx.customModelProperties.register(NtCTubeProvider, this.params.autoAttach);

            this.ctx.representation.structure.themes.colorThemeRegistry.add(ConfalPyramidsColorThemeProvider);
            this.ctx.representation.structure.registry.add(ConfalPyramidsRepresentationProvider);
            this.ctx.representation.structure.themes.colorThemeRegistry.add(NtCTubeColorThemeProvider);
            this.ctx.representation.structure.registry.add(NtCTubeRepresentationProvider);

            this.ctx.builders.structure.representation.registerPreset(ConfalPyramidsPreset);
            this.ctx.builders.structure.representation.registerPreset(NtCTubePreset);
        }

        unregister() {
            this.ctx.customModelProperties.unregister(ConfalPyramidsProvider.descriptor.name);
            this.ctx.customModelProperties.unregister(NtCTubeProvider.descriptor.name);

            this.ctx.representation.structure.registry.remove(ConfalPyramidsRepresentationProvider);
            this.ctx.representation.structure.themes.colorThemeRegistry.remove(ConfalPyramidsColorThemeProvider);
            this.ctx.representation.structure.registry.remove(NtCTubeRepresentationProvider);
            this.ctx.representation.structure.themes.colorThemeRegistry.remove(NtCTubeColorThemeProvider);

            this.ctx.builders.structure.representation.unregisterPreset(ConfalPyramidsPreset);
            this.ctx.builders.structure.representation.unregisterPreset(NtCTubePreset);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(true),
        showToolTip: PD.Boolean(true)
    })
});

