/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ConfalPyramidsColorThemeProvider } from './color';
import { ConfalPyramids, ConfalPyramidsProvider } from './property';
import { ConfalPyramidsRepresentationProvider } from './representation';
import { ConfalPyramidsTypes } from './types';
import { PluginBehavior } from '../../../mol-plugin/behavior/behavior';
import { StructureRepresentationPresetProvider, PresetStructureRepresentations } from '../../../mol-plugin-state/builder/structure/representation-preset';
import { StateObjectRef } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

export const DnatcoConfalPyramidsPreset = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-confal-pyramids',
    display: {
        name: 'Confal Pyramids', group: 'Annotation',
        description: 'Schematic depiction of conformer class and confal value.',
    },
    isApplicable(a) {
        return a.data.models.length >= 1 && a.data.models.some(m => ConfalPyramids.isApplicable(m));
    },
    params: () => StructureRepresentationPresetProvider.CommonParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        const model = structureCell?.obj?.data.model;
        if (!structureCell || !model) return {};

        await plugin.runTask(Task.create('Confal Pyramids', async runtime => {
            await ConfalPyramidsProvider.attach({ runtime, assetManager: plugin.managers.asset }, model);
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

export const DnatcoConfalPyramids = PluginBehavior.create<{ autoAttach: boolean, showToolTip: boolean }>({
    name: 'dnatco-confal-pyramids-prop',
    category: 'custom-props',
    display: {
        name: 'Confal Pyramids',
        description: 'Schematic depiction of conformer class and confal value.',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showToolTip: boolean }> {
        private provider = ConfalPyramidsProvider;

        register(): void {
            this.ctx.customModelProperties.register(this.provider, this.params.autoAttach);

            this.ctx.representation.structure.themes.colorThemeRegistry.add(ConfalPyramidsColorThemeProvider);
            this.ctx.representation.structure.registry.add(ConfalPyramidsRepresentationProvider);

            this.ctx.builders.structure.representation.registerPreset(DnatcoConfalPyramidsPreset);
        }

        update(p: { autoAttach: boolean, showToolTip: boolean }) {
            const updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.params.showToolTip = p.showToolTip;
            this.ctx.customModelProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        unregister() {
            this.ctx.customModelProperties.unregister(ConfalPyramidsProvider.descriptor.name);

            this.ctx.representation.structure.registry.remove(ConfalPyramidsRepresentationProvider);
            this.ctx.representation.structure.themes.colorThemeRegistry.remove(ConfalPyramidsColorThemeProvider);

            this.ctx.builders.structure.representation.unregisterPreset(DnatcoConfalPyramidsPreset);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(true),
        showToolTip: PD.Boolean(true)
    })
});

export function confalPyramidLabel(halfPyramid: ConfalPyramidsTypes.HalfPyramid) {
    const { step } = halfPyramid;
    return `
        <b>${step.auth_asym_id_1}</b> |
        <b>${step.label_comp_id_1} ${step.auth_seq_id_1}${step.PDB_ins_code_1}${step.label_alt_id_1.length > 0 ? ` (alt ${step.label_alt_id_1})` : ''}
           ${step.label_comp_id_2} ${step.auth_seq_id_2}${step.PDB_ins_code_2}${step.label_alt_id_2.length > 0 ? ` (alt ${step.label_alt_id_2})` : ''} </b><br />
        <i>NtC:</i> ${step.NtC} | <i>Confal score:</i> ${step.confal_score} | <i>RMSD:</i> ${step.rmsd.toFixed(2)}
    `;
}
