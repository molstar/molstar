import {
    PresetStructureRepresentations,
    StructureRepresentationPresetProvider,
} from '../../../mol-plugin-state/builder/structure/representation-preset';
import { StateObjectRef } from '../../../mol-state';
import { SbNcbrPartialChargesPropertyProvider } from './property';
import { SbNcbrPartialChargesColorThemeProvider } from './color';

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
        return PresetStructureRepresentations.auto.apply(
            ref,
            { ...params, theme: { globalName: colorTheme, focus: { name: colorTheme, params: { chargeType: 'atom' } } } },
            plugin
        );
    },
});
