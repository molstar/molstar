/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ConfalPyramidsPreset } from '@molstar/dnatco-extension/confal-pyramids/behavior';
import { ConfalPyramidsColorThemeProvider } from '@molstar/dnatco-extension/confal-pyramids/color';
import { ConfalPyramidsProvider } from '@molstar/dnatco-extension/confal-pyramids/property';
import { ConfalPyramidsRepresentationProvider } from '@molstar/dnatco-extension/confal-pyramids/representation';
import { NtCTubePreset } from '@molstar/dnatco-extension/ntc-tube/behavior';
import { NtCTubeColorThemeProvider } from '@molstar/dnatco-extension/ntc-tube/color';
import { NtCTubeProvider } from '@molstar/dnatco-extension/ntc-tube/property';
import { NtCTubeRepresentationProvider } from '@molstar/dnatco-extension/ntc-tube/representation';

export const DnatcoNtCs = PluginBehavior.create<{ autoAttach: boolean; showToolTip: boolean }>({
  name: 'dnatco-ntcs',
  category: 'custom-props',
  display: {
    name: 'DNATCO NtC Annotations',
    description: 'DNATCO NtC Annotations',
  },
  ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean; showToolTip: boolean }> {
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
    showToolTip: PD.Boolean(true),
  }),
});
