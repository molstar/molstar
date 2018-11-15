/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateAction } from 'mol-state/action';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { StateSelection } from 'mol-state/state/selection';
import { CartoonParams } from 'mol-repr/structure/representation/cartoon';

export const CreateStructureFromPDBe = StateAction.create<PluginStateObject.Root, void, { id: string }>({
    from: [PluginStateObject.Root],
    display: {
        name: 'Entry from PDBe',
        description: 'Download a structure from PDBe and create its default Assembly and visual'
    },
    params: {
        default: () => ({ id: '1grm' }),
        definition: () => ({
            id: PD.Text('1grm', { label: 'PDB id' }),
        }),
        // validate: p => !p.id || !p.id.trim() ? [['Enter id.', 'id']] : void 0
    },
    apply({ params, state }) {
        const url = `http://www.ebi.ac.uk/pdbe/static/entry/${params.id.toLowerCase()}_updated.cif`;
        const b = state.build();

        // const query = MolScriptBuilder.struct.generator.atomGroups({
        //     // 'atom-test': MolScriptBuilder.core.rel.eq([
        //     //     MolScriptBuilder.struct.atomProperty.macromolecular.label_comp_id(),
        //     //     MolScriptBuilder.es('C')
        //     // ]),
        //     'residue-test': MolScriptBuilder.core.rel.eq([
        //         MolScriptBuilder.struct.atomProperty.macromolecular.label_comp_id(),
        //         'ALA'
        //     ])
        // });

        const newTree = b.toRoot()
            .apply(StateTransforms.Data.Download, { url })
            .apply(StateTransforms.Data.ParseCif)
            .apply(StateTransforms.Model.ParseTrajectoryFromMmCif, {})
            .apply(StateTransforms.Model.CreateModelFromTrajectory, { modelIndex: 0 })
            .apply(StateTransforms.Model.CreateStructureAssembly)
            // .apply(StateTransforms.Model.CreateStructureSelection, { query, label: 'ALA residues' })
            .apply(StateTransforms.Visuals.CreateStructureRepresentation, {
                type: {
                    name: 'cartoon',
                    params: PD.getDefaultValues(CartoonParams)
                }
            })
            .getTree();

        return state.update(newTree);
    }
});

export const UpdateTrajectory = StateAction.create<PluginStateObject.Root, void, { action: 'advance' | 'reset', by?: number }>({
    from: [],
    display: {
        name: 'Update Trajectory'
    },
    params: {
        default: () => ({ action: 'reset', by: 1 })
    },
    apply({ params, state }) {
        const models = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Model).filter(c => c.transform.transformer === StateTransforms.Model.CreateModelFromTrajectory));

        const update = state.build();

        if (params.action === 'reset') {
            for (const m of models) {
                update.to(m.transform.ref).update(StateTransforms.Model.CreateModelFromTrajectory,
                    () => ({ modelIndex: 0}));
            }
        } else {
            for (const m of models) {
                const parent = StateSelection.findAncestorOfType(state.tree, state.cells, m.transform.ref, [PluginStateObject.Molecule.Trajectory]);
                if (!parent || !parent.obj) continue;
                const traj = parent.obj as PluginStateObject.Molecule.Trajectory;
                update.to(m.transform.ref).update(StateTransforms.Model.CreateModelFromTrajectory,
                    old => {
                        let modelIndex = (old.modelIndex + params.by!) % traj.data.length;
                        if (modelIndex < 0) modelIndex += traj.data.length;
                        return { modelIndex };
                    });
            }
        }

        return state.update(update);
    }
});