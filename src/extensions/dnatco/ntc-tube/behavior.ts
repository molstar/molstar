/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { NtCTubeColorThemeProvider } from './color';
import { NtCTubeProvider } from './property';
import { NtCTubeRepresentationProvider } from './representation';
import { DnatcoTypes } from '../types';
import { Dnatco } from '../property';
import { StructureRepresentationPresetProvider, PresetStructureRepresentations } from '../../../mol-plugin-state/builder/structure/representation-preset';
import { StateObjectRef } from '../../../mol-state';
import { Task } from '../../../mol-task';

export const NtCTubePreset = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-ntc-tube',
    display: {
        name: 'NtC Tube', group: 'Annotation',
        description: 'NtC Tube',
    },
    isApplicable(a) {
        return a.data.models.length >= 1 && a.data.models.some(m => Dnatco.isApplicable(m));
    },
    params: () => StructureRepresentationPresetProvider.CommonParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        const model = structureCell?.obj?.data.model;
        if (!structureCell || !model) return {};

        await plugin.runTask(Task.create('NtC tube', async runtime => {
            await NtCTubeProvider.attach({ runtime, assetManager: plugin.managers.asset, errorContext: plugin.errorContext }, model);
        }));

        const { components, representations } = await PresetStructureRepresentations.auto.apply(ref, { ...params }, plugin);

        const tube = await plugin.builders.structure.tryCreateComponentStatic(structureCell, 'nucleic', { label: 'NtC Tube' });
        const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, params);

        let tubeRepr;
        if (representations)
            tubeRepr = builder.buildRepresentation(update, tube, { type: NtCTubeRepresentationProvider, typeParams, color: NtCTubeColorThemeProvider }, { tag: 'ntc-tube' });

        await update.commit({ revertOnError: true });
        return { components: { ...components, tube }, representations: { ...representations, tubeRepr } };
    }
});

export function NtCTubeSegmentLabel(step: DnatcoTypes.Step) {
    return `
        <b>${step.auth_asym_id_1}</b> |
        <b>${step.label_comp_id_1} ${step.auth_seq_id_1}${step.PDB_ins_code_1}${step.label_alt_id_1.length > 0 ? ` (alt ${step.label_alt_id_1})` : ''}
           ${step.label_comp_id_2} ${step.auth_seq_id_2}${step.PDB_ins_code_2}${step.label_alt_id_2.length > 0 ? ` (alt ${step.label_alt_id_2})` : ''} </b><br />
        <i>NtC:</i> ${step.NtC} | <i>Confal score:</i> ${step.confal_score} | <i>RMSD:</i> ${step.rmsd.toFixed(2)}
    `;
}
