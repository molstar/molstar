/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransforms } from '../../transforms';
import { StructureComplexElementTypes } from '../../transforms/model';
import { StructureRepresentation3DHelpers } from '../../transforms/representation';
import { StructureSelectionQueries as Q } from '../../../util/structure-selection-helper';
import { BuiltInStructureRepresentations } from '../../../../mol-repr/structure/registry';
import { StructureRepresentationProvider, RepresentationProviderTags } from './provider';
import { StateBuilder } from '../../../../mol-state';
import { PluginStateObject } from '../../objects';

const auto = StructureRepresentationProvider({
    id: 'preset-structure-representation-auto',
    display: { name: 'Automaic', group: 'Preset' },
    apply(state, structureCell, _, plugin) {
        const s = structureCell.obj!.data;

        // TODO: a way to improve this?

        if (s.elementCount < 50000) {
            return defaultPreset.apply(state, structureCell, void 0, plugin);
        } else if (s.elementCount < 200000) {
            return proteinAndNucleic.apply(state, structureCell, void 0, plugin);
        } else {
            if (s.unitSymmetryGroups[0].units.length > 10) {
                return capsid.apply(state, structureCell, void 0, plugin);
            } else {
                return coarseCapsid.apply(state, structureCell, void 0, plugin);
            }
        }
    }
});

const defaultPreset = StructureRepresentationProvider({
    id: 'preset-structure-representation-default',
    display: { name: 'Default', group: 'Preset' },
    apply(state, structureCell, _, plugin) {
        const root = state.build().to(structureCell.transform.ref);
        const structure = structureCell.obj!.data;

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        applyComplex(root, 'protein-or-nucleic')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParams(plugin, 'cartoon', structure));

        applyComplex(root, 'ligand')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure));

        applyComplex(root, 'modified')
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

        return state.updateTree(root, { revertIfAborted: true });
    }
});

const proteinAndNucleic = StructureRepresentationProvider({
    id: 'preset-structure-representation-protein-and-nucleic',
    display: { name: 'Protein & Nucleic', group: 'Preset' },
    apply(state, structureCell, _, plugin) {
        const root = plugin.state.dataState.build().to(structureCell.transform.ref);
        const structure = structureCell.obj!.data;
        const reprTags = [this.id, RepresentationProviderTags.Representation];

        applySelection(root, 'protein')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParams(plugin, 'cartoon', structure));

        applySelection(root, 'nucleic')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParams(plugin, 'gaussian-surface', structure));

        return state.updateTree(root, { revertIfAborted: true });
    }
});

const capsid = StructureRepresentationProvider({
    id: 'preset-structure-representation-capsid',
    display: { name: 'Capsid', group: 'Preset' },
    apply(state, structureCell, _, plugin) {
        const root = plugin.state.dataState.build().to(structureCell.transform.ref);
        const structure = structureCell.obj!.data;

        const params = StructureRepresentation3DHelpers.createParams(plugin, structure, {
            repr: [BuiltInStructureRepresentations['gaussian-surface'], () => ({ smoothness: 1 })]
        });

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        applySelection(root, 'polymer')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, params);

        return state.updateTree(root, { revertIfAborted: true });
    }
});

const coarseCapsid = StructureRepresentationProvider({
    id: 'preset-structure-representation-coarse-capsid',
    display: { name: 'Coarse Capsid', group: 'Preset' },
    apply(state, structureCell, _, plugin) {
        const root = plugin.state.dataState.build().to(structureCell.transform.ref);
        const structure = structureCell.obj!.data;

        const params = StructureRepresentation3DHelpers.createParams(plugin, structure, {
            repr: [
                BuiltInStructureRepresentations['gaussian-surface'],
                () => ({ smoothness: 0.5, radiusOffset: 1, /*visuals: ['gaussian-surface-mesh']*/ })
            ]
        });

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        applySelection(root, 'trace')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, params);

        return state.updateTree(root, { revertIfAborted: true });
    }
});

const cartoon = StructureRepresentationProvider({
    id: 'preset-structure-representation-cartoon',
    display: { name: 'Cartoon', group: 'Preset' },
    apply(state, structureCell, _, plugin) {
        const root = plugin.state.dataState.build().to(structureCell.transform.ref);
        const structure = structureCell.obj!.data;

        const params = StructureRepresentation3DHelpers.createParams(plugin, structure, {
            repr: BuiltInStructureRepresentations['cartoon']
        });

        const reprTags = [this.id, RepresentationProviderTags.Representation];

        applySelection(root, 'polymer')
            .applyOrUpdateTagged(reprTags, StateTransforms.Representation.StructureRepresentation3D, params);

        return state.updateTree(root, { revertIfAborted: true });
    }
});

function applyComplex(to: StateBuilder.To<PluginStateObject.Molecule.Structure>, type: keyof typeof StructureComplexElementTypes) {
    return to.applyOrUpdateTagged(type, StateTransforms.Model.StructureComplexElement, { type }, { tags: RepresentationProviderTags.Selection });
}

function applySelection(to: StateBuilder.To<PluginStateObject.Molecule.Structure>, query: keyof typeof Q) {
    return to.applyOrUpdateTagged(query, StateTransforms.Model.StructureSelectionFromExpression,
        { expression: Q[query].expression, label: Q[query].label },
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