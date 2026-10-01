/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ConfalPyramidsColorThemeProvider } from './color.js';
import { ConfalPyramidsProvider } from './property.js';
import { ConfalPyramidsRepresentationProvider } from './representation.js';
import { Dnatco } from '@molstar/dnatco-extension/property';
import { DnatcoTypes } from '@molstar/dnatco-extension/types';
import { StructureRepresentationPresetProvider, PresetStructureRepresentations } from '@molstar/plugin/state/builder/structure/representation-preset';
import { StateObjectRef } from '@molstar/core/state';
import { Task } from '@molstar/core/task';

export const ConfalPyramidsPreset = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-confal-pyramids',
    display: {
        name: 'Confal Pyramids', group: 'Annotation',
        description: 'Schematic depiction of conformer class and confal value.',
    },
    isApplicable(a) {
        return a.data.models.length >= 1 && a.data.models.some(m => Dnatco.isApplicable(m));
    },
    params: () => StructureRepresentationPresetProvider.CommonParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        const model = structureCell?.obj?.data.model;
        if (!structureCell || !model) return {};

        await plugin.runTask(Task.create('Confal Pyramids', async runtime => {
            await ConfalPyramidsProvider.attach({ runtime, assetManager: plugin.managers.asset, errorContext: plugin.errorContext }, model);
        }));

        const { components, representations } = await PresetStructureRepresentations.auto.apply(ref, { ...params }, plugin);

        const pyramids = await plugin.builders.structure.tryCreateComponentStatic(structureCell, 'nucleic', { label: 'Confal Pyramids' });
        const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, params);

        let pyramidsRepr;
        if (representations)
            pyramidsRepr = builder.buildRepresentation(update, pyramids, { type: ConfalPyramidsRepresentationProvider, typeParams, color: ConfalPyramidsColorThemeProvider }, { tag: 'confal-pyramdis' });

        await update.commit({ revertOnError: true });
        return { components: { ...components, pyramids }, representations: { ...representations, pyramidsRepr } };
    }
});

const RemoveNewline = /\r?\n/g;
export function confalPyramidLabel(step: DnatcoTypes.Step) {
    return `
        <b>${step.auth_asym_id_1}</b> |
        <b>${step.label_comp_id_1} ${step.auth_seq_id_1}${step.PDB_ins_code_1}${step.label_alt_id_1.length > 0 ? ` (alt ${step.label_alt_id_1})` : ''}
           ${step.label_comp_id_2} ${step.auth_seq_id_2}${step.PDB_ins_code_2}${step.label_alt_id_2.length > 0 ? ` (alt ${step.label_alt_id_2})` : ''} </b><br />
        <i>NtC:</i> ${step.NtC} | <i>Confal score:</i> ${step.confal_score} | <i>RMSD:</i> ${step.rmsd.toFixed(2)}
    `.replace(RemoveNewline, '');
}
