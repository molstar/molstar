/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransforms } from '../../transforms';
import { StructureRepresentation3DHelpers } from '../../transforms/representation';
import { StructureSelectionQueries as Q } from '../../../mol-plugin/util/structure-selection-helper';
import { BuiltInStructureRepresentations } from '../../../mol-repr/structure/registry';
import { StructureRepresentationProvider, RepresentationProviderTags } from './provider';
import { StateObjectRef } from '../../../mol-state';
import { PluginStateObject } from '../../objects';
import { StaticStructureComponentType } from '../../helpers/structure-component';
import { PluginContext } from '../../../mol-plugin/context';

const auto = StructureRepresentationProvider({
    id: 'preset-structure-representation-auto',
    display: { name: 'Automaic', group: 'Preset' },
    apply(ctx, state, structureCell, _, plugin) {
        const s = structureCell.obj!.data;

        // TODO: a way to improve this?

        if (s.elementCount < 50000) {
            return defaultPreset.apply(ctx, state, structureCell, void 0, plugin);
        } else if (s.elementCount < 200000) {
            return proteinAndNucleic.apply(ctx, state, structureCell, void 0, plugin);
        } else {
            if (s.unitSymmetryGroups[0].units.length > 10) {
                return capsid.apply(ctx, state, structureCell, void 0, plugin);
            } else {
                return coarseCapsid.apply(ctx, state, structureCell, void 0, plugin);
            }
        }
    }
});

const defaultPreset = StructureRepresentationProvider({
    id: 'preset-structure-representation-default',
    display: { name: 'Default', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const structure = structureCell.obj?.data!;
        const reprTags = [this.id, RepresentationProviderTags.Representation];
        
        const components = {
            proteinOrNucleic: await staticComponent(plugin, structureCell, 'protein-or-nucleic'),
            ligand: await staticComponent(plugin, structureCell, 'ligand'),
            nonStandard: await staticComponent(plugin, structureCell, 'non-standard'),
            branched: await staticComponent(plugin, structureCell, 'branched'),
            water: await staticComponent(plugin, structureCell, 'water'),
            coarse: await staticComponent(plugin, structureCell, 'coarse')
        };
                
        const builder = state.build();
        const representations = {
            proteinOrNucleic: components.proteinOrNucleic && builder
                .to(components.proteinOrNucleic)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'cartoon', structure)).selector,
            ligand: components.ligand && builder
                .to(components.ligand)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure)).selector,
            nonStandard: components.nonStandard && builder
                .to(components.nonStandard)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParamsWithTheme(plugin, 'ball-and-stick', 'polymer-id', structure, void 0)).selector,
            branched: components.branched && {
                ballAndStick: builder.to(components.branched).applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure, { alpha: 0.15 })).selector,
                snfg3d: builder.to(components.branched).applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'carbohydrate', structure)).selector
            },
            water: components.water && builder
                .to(components.water)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure, { alpha: 0.51 })).selector,
            coarse: components.coarse && builder
                .to(components.coarse)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParamsWithTheme(plugin, 'spacefill', 'polymer-id', structure, {}))
        };
        
        await state.updateTree(builder, { revertOnError: false }).runInContext(ctx);
        return { components, representations };
    }
});

const proteinAndNucleic = StructureRepresentationProvider({
    id: 'preset-structure-representation-protein-and-nucleic',
    display: { name: 'Protein & Nucleic', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const structure = structureCell.obj!.data;
        const reprTags = [this.id, RepresentationProviderTags.Representation];

        const components = {
            protein: await selectionComponent(plugin, structureCell, 'protein'),
            nucleic: await selectionComponent(plugin, structureCell, 'nucleic'),
        };

        const builder = state.build();
        const representations = {
            protein: components.protein && builder
                .to(components.protein)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'cartoon', structure)).selector,
            nucleic: components.nucleic && builder
                .to(components.nucleic)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, StructureRepresentation3DHelpers.getDefaultParams(plugin, 'gaussian-surface', structure)).selector,
        };
      
        await state.updateTree(builder, { revertOnError: true }).runInContext(ctx);
        return { components, representations };
    }
});

const capsid = StructureRepresentationProvider({
    id: 'preset-structure-representation-capsid',
    display: { name: 'Capsid', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const structure = structureCell.obj!.data;

        const params = StructureRepresentation3DHelpers.createParams(plugin, structure, {
            repr: [BuiltInStructureRepresentations['gaussian-surface'], () => ({ smoothness: 1 })]
        });

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        const components = {
            polymer: await selectionComponent(plugin, structureCell, 'polymer')
        };

        const builder = state.build();
        const representations = {
            polymer: components.polymer && builder
                .to(components.polymer)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, params).selector,
        };

        await state.updateTree(builder, { revertOnError: true }).runInContext(ctx);
        return { components, representations };
    }
});

const coarseCapsid = StructureRepresentationProvider({
    id: 'preset-structure-representation-coarse-capsid',
    display: { name: 'Coarse Capsid', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const structure = structureCell.obj!.data;

        const params = StructureRepresentation3DHelpers.createParams(plugin, structure, {
            repr: [
                BuiltInStructureRepresentations['gaussian-surface'],
                () => ({ smoothness: 0.5, radiusOffset: 1, /* visuals: ['gaussian-surface-mesh']*/ })
            ]
        });

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        const components = {
            trace: await selectionComponent(plugin, structureCell, 'trace')
        };

        const builder = state.build();
        const representations = {
            trace: components.trace && builder
                .to(components.trace)
                .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, params).selector,
        };
        
        await state.updateTree(builder, { revertOnError: true }).runInContext(ctx);

        return { components, representations };
    }
});

function staticComponent(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, type: StaticStructureComponentType) {
    return plugin.builders.structure.tryCreateComponent(structure, { 
        type: { name: 'static', params: type },
        nullIfEmpty: true,
        label: ''
    }, `static-${type}`);
}

function selectionComponent(plugin: PluginContext, structure: StateObjectRef<PluginStateObject.Molecule.Structure>, query: keyof typeof Q) {
    return plugin.builders.structure.tryCreateComponent(structure, { 
        type: { name: 'expression', params: Q[query].expression },
        nullIfEmpty: true,
        label: Q[query].label
    }, `selection-${query}`);
}

export const PresetStructureReprentations = {
    auto,
    default: defaultPreset,
    proteinAndNucleic,
    capsid,
    coarseCapsid
};
export type PresetStructureReprentations = typeof PresetStructureReprentations;