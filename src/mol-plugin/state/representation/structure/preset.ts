/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransforms } from '../../transforms';
import { StructureComplexElementTypes } from '../../transforms/model';
import { StructureRepresentation3DHelpers } from '../../transforms/representation';
import { applyBuiltInSelection } from '../../../util/structure-selection-helper';
import { BuiltInStructureRepresentations } from '../../../../mol-repr/structure/registry';
import { StructureRepresentationProvider } from './providers';

export const PresetStructureReprentations = {
    default: StructureRepresentationProvider({
        id: 'preset-structure-representation-default',
        display: { name: 'Default', group: 'Preset' },
        apply(state, structureCell, _, plugin) {
            const root = state.build().to(structureCell.transform.ref);
            const structure = structureCell.obj!.data;

            const tags = this.id;

            root.apply(StateTransforms.Model.StructureComplexElement, { type: 'protein-or-nucleic' }, { tags: [this.id, StructureComplexElementTypes['protein-or-nucleic']] })
                .apply(StateTransforms.Representation.StructureRepresentation3D,
                    StructureRepresentation3DHelpers.getDefaultParams(plugin, 'cartoon', structure), { tags });

            root.apply(StateTransforms.Model.StructureComplexElement, { type: 'ligand' }, { tags: [this.id, StructureComplexElementTypes.ligand] })
                .apply(StateTransforms.Representation.StructureRepresentation3D,
                    StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure), { tags });

            root.apply(StateTransforms.Model.StructureComplexElement, { type: 'modified' }, { tags: [this.id, StructureComplexElementTypes.modified] })
                .apply(StateTransforms.Representation.StructureRepresentation3D,
                    StructureRepresentation3DHelpers.getDefaultParamsWithTheme(plugin, 'ball-and-stick', 'polymer-id', structure, void 0), { tags });

            const branched = root.apply(StateTransforms.Model.StructureComplexElement, { type: 'branched' }, { tags: [this.id, StructureComplexElementTypes.branched] })

            branched.apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure, { alpha: 0.15 }), { tags });
            branched.apply(StateTransforms.Representation.StructureRepresentation3D,
                StructureRepresentation3DHelpers.getDefaultParams(plugin, 'carbohydrate', structure), { tags });

            root.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' }, { tags: [this.id, StructureComplexElementTypes.water] })
                .apply(StateTransforms.Representation.StructureRepresentation3D,
                    StructureRepresentation3DHelpers.getDefaultParams(plugin, 'ball-and-stick', structure, { alpha: 0.51 }), { tags });

            root.apply(StateTransforms.Model.StructureComplexElement, { type: 'coarse' }, { tags: [this.id, StructureComplexElementTypes.coarse] })
                .apply(StateTransforms.Representation.StructureRepresentation3D,
                    StructureRepresentation3DHelpers.getDefaultParamsWithTheme(plugin, 'spacefill', 'polymer-id', structure, {}), { tags });

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

            applyBuiltInSelection(root, 'polymer', this.id)
                .apply(StateTransforms.Representation.StructureRepresentation3D, params, { tags: this.id });

            return state.updateTree(root, { revertIfAborted: true });
        }
    })
};