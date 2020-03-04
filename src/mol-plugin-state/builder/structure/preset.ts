/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransforms } from '../../transforms';
import { StructureRepresentation3DHelpers } from '../../transforms/representation';
import { StructureSelectionQueries as Q } from '../../../mol-plugin/util/structure-selection-helper';
import { BuiltInStructureRepresentations } from '../../../mol-repr/structure/registry';
import { StructureRepresentationProvider, RepresentationProviderTags } from './provider';
import { StateBuilder } from '../../../mol-state';
import { PluginStateObject } from '../../objects';
import { StaticStructureComponentType } from '../../helpers/structure-component';

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
        const root = state.build().to(structureCell.transform.ref);
        const structure = structureCell.obj!.data;

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        applyComplex(root, 'protein-or-nucleic')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParams(plugin, 'cartoon', structure));

        const ligand = applyComplex(root, 'ligand');
        const ligandRepr = ligand.applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure));

        applyComplex(root, 'non-standard')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsWithTheme(plugin, 'ball-and-stick', 'polymer-id', structure, void 0));

        const branched = applyComplex(root, 'branched');

        branched.applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure, { alpha: 0.15 }));
        branched.applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParams(plugin, 'carbohydrate', structure));

        applyComplex(root, 'water')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure, { alpha: 0.51 }));

        applyComplex(root, 'coarse')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParamsWithTheme(plugin, 'spacefill', 'polymer-id', structure, {}));

        await state.updateTree(root, { revertOnError: true }).runInContext(ctx);

        return {
            ligand: {
                selection: ligand.selector,
                repr: ligandRepr.selector
            }
        };
    }
});

const proteinAndNucleic = StructureRepresentationProvider({
    id: 'preset-structure-representation-protein-and-nucleic',
    display: { name: 'Protein & Nucleic', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const root = plugin.state.dataState.build().to(structureCell.transform.ref);
        const structure = structureCell.obj!.data;
        const reprTags = [this.id, RepresentationProviderTags.Representation];

        applySelection(root, 'protein')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParams(plugin, 'cartoon', structure));

        applySelection(root, 'nucleic')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParams(plugin, 'gaussian-surface', structure));

        await state.updateTree(root, { revertOnError: true }).runInContext(ctx);
        return {};
    }
});

const capsid = StructureRepresentationProvider({
    id: 'preset-structure-representation-capsid',
    display: { name: 'Capsid', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const root = plugin.state.dataState.build().to(structureCell.transform.ref);
        const structure = structureCell.obj!.data;

        const params = StructureRepresentation3DHelpers.createParams(plugin, structure, {
            repr: [BuiltInStructureRepresentations['gaussian-surface'], () => ({ smoothness: 1 })]
        });

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        applySelection(root, 'polymer')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, params);

        await state.updateTree(root, { revertOnError: true }).runInContext(ctx);
        return {};
    }
});

const coarseCapsid = StructureRepresentationProvider({
    id: 'preset-structure-representation-coarse-capsid',
    display: { name: 'Coarse Capsid', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const root = plugin.state.dataState.build().to(structureCell.transform.ref);
        const structure = structureCell.obj!.data;

        const params = StructureRepresentation3DHelpers.createParams(plugin, structure, {
            repr: [
                BuiltInStructureRepresentations['gaussian-surface'],
                () => ({ smoothness: 0.5, radiusOffset: 1, /* visuals: ['gaussian-surface-mesh']*/ })
            ]
        });

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        applySelection(root, 'trace')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, params);

        await state.updateTree(root, { revertOnError: true }).runInContext(ctx);
        return {};
    }
});

const cartoon = StructureRepresentationProvider({
    id: 'preset-structure-representation-cartoon',
    display: { name: 'Cartoon', group: 'Preset' },
    async apply(ctx, state, structureCell, _, plugin) {
        const root = plugin.state.dataState.build().to(structureCell.transform.ref);
        const structure = structureCell.obj!.data;

        const params = StructureRepresentation3DHelpers.createParams(plugin, structure, {
            repr: BuiltInStructureRepresentations['cartoon']
        });

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        applySelection(root, 'polymer')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, params);

        await state.updateTree(root, { revertOnError: true }).runInContext(ctx);
        return {};
    }
});

function applyComplex(to: StateBuilder.To<PluginStateObject.Molecule.Structure>, type: StaticStructureComponentType) {
    return to.applyOrUpdateTagged(type, StateTransforms.Model.StructureComponent, { 
        type: { name: 'static', params: type },
        nullIfEmpty: true,
        label: ''
    }, { tags: RepresentationProviderTags.Selection });
}

function applySelection(to: StateBuilder.To<PluginStateObject.Molecule.Structure>, query: keyof typeof Q) {
    return to.applyOrUpdateTagged(query, StateTransforms.Model.StructureComponent, { 
        type: { name: 'expression', params: Q[query].expression },
        nullIfEmpty: true,
        label: Q[query].label
    },
    { tags: RepresentationProviderTags.Selection });
}

export const PresetStructureReprentations = {
    auto,
    default: defaultPreset,
    proteinAndNucleic,
    capsid,
    coarseCapsid,
    cartoon
};
export type PresetStructureReprentations = typeof PresetStructureReprentations;