import type { LociLabelProvider } from '@molstar/plugin/state/manager/loci-label';
import { PluginBehavior } from '@molstar/plugin/behavior';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { SbNcbrPartialChargesColorThemeProvider } from './color.js';
import { SbNcbrPartialChargesPropertyProvider } from './property.js';
import { SbNcbrPartialChargesLociLabelProvider } from './labels.js';
import { SbNcbrPartialChargesPreset } from './preset.js';

export const SbNcbrPartialCharges = PluginBehavior.create<{ autoAttach: boolean; showToolTip: boolean }>({
  name: 'sb-ncbr-partial-charges',
  category: 'misc',
  display: {
    name: 'SB NCBR Partial Charges',
  },
  ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean; showToolTip: boolean }> {
    private SbNcbrPartialChargesLociLabelProvider: LociLabelProvider = SbNcbrPartialChargesLociLabelProvider(this.ctx);

    register(): void {
      this.ctx.customModelProperties.register(SbNcbrPartialChargesPropertyProvider, this.params.autoAttach);
      this.ctx.representation.structure.themes.colorThemeRegistry.add(SbNcbrPartialChargesColorThemeProvider);
      this.ctx.managers.lociLabels.addProvider(this.SbNcbrPartialChargesLociLabelProvider);
      this.ctx.builders.structure.representation.registerPreset(SbNcbrPartialChargesPreset);
    }

    unregister() {
      this.ctx.customModelProperties.unregister(SbNcbrPartialChargesPropertyProvider.descriptor.name);
      this.ctx.representation.structure.themes.colorThemeRegistry.remove(SbNcbrPartialChargesColorThemeProvider);
      this.ctx.managers.lociLabels.removeProvider(this.SbNcbrPartialChargesLociLabelProvider);
      this.ctx.builders.structure.representation.unregisterPreset(SbNcbrPartialChargesPreset);
    }
  },
  params: () => ({
    autoAttach: PD.Boolean(true),
    showToolTip: PD.Boolean(true),
  }),
});
