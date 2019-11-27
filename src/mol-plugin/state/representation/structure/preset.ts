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

export const PresetStructureReprentations = {
    default: StructureRepresentationProvider({
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
    }),
    capsid: StructureRepresentationProvider({
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
    }),
    cartoon: StructureRepresentationProvider({
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
    })
};

function applyComplex(to: StateBuilder.To<PluginStateObject.Molecule.Structure>, type: keyof typeof StructureComplexElementTypes) {
    return to.applyOrUpdateTagged(type, StateTransforms.Model.StructureComplexElement, { type }, { tags: RepresentationProviderTags.Selection });
}

function applySelection(to: StateBuilder.To<PluginStateObject.Molecule.Structure>, query: keyof typeof Q) {
    return to.applyOrUpdateTagged(query, StateTransforms.Model.StructureSelectionFromExpression,
        { expression: Q[query].expression, label: Q[query].label },
        { tags: RepresentationProviderTags.Selection });
}