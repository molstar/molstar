import type { LociLabelProvider } from '@molstar/plugin/state/manager/loci-label';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
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

    private unregisterEntry: (() => void) | undefined;

    register(): void {
      this.ctx.customModelProperties.register(SbNcbrPartialChargesPropertyProvider, this.params.autoAttach);

      const entry: PluginRegistryEntry = {
        structure: {
          themes: { color: [SbNcbrPartialChargesColorThemeProvider] },
          presets: { representation: [SbNcbrPartialChargesPreset] },
        },
        lociLabels: [this.SbNcbrPartialChargesLociLabelProvider],
      };
      this.unregisterEntry = this.ctx.register(entry);
    }

    unregister() {
      this.ctx.customModelProperties.unregister(SbNcbrPartialChargesPropertyProvider.descriptor.name);
      this.unregisterEntry?.();
      this.unregisterEntry = undefined;
    }
  },
  params: () => ({
    autoAttach: PD.Boolean(true),
    showToolTip: PD.Boolean(true),
  }),
});
