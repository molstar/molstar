import { StructureRepresentationPresetProvider } from '@molstar/plugin/state/builder/structure/representation-presets/types';
import { AutoPreset } from '@molstar/plugin/state/builder/structure/representation-presets/auto';
import { StateObjectRef } from '@molstar/core/state';
import { SbNcbrPartialChargesPropertyProvider } from './property.js';
import { SbNcbrPartialChargesColorThemeProvider } from './color.js';

export const SbNcbrPartialChargesPreset = StructureRepresentationPresetProvider({
  id: 'sb-ncbr-partial-charges-preset',
  display: {
    name: 'SB NCBR Partial Charges',
    group: 'Annotation',
    description: 'Color atoms and residues based on their partial charge.',
  },
  isApplicable(a) {
    return !!a.data.models.some((m) => SbNcbrPartialChargesPropertyProvider.isApplicable(m));
  },
  params: () => StructureRepresentationPresetProvider.CommonParams,
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    const structure = structureCell?.obj?.data;
    if (!structureCell || !structure) return {};

    const colorTheme = SbNcbrPartialChargesColorThemeProvider.name as any;
    return AutoPreset.apply(
      ref,
      { ...params, theme: { globalName: colorTheme, focus: { name: colorTheme, params: { chargeType: 'atom' } } } },
      plugin,
    );
  },
});
