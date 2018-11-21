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
import { BallAndStickParams } from 'mol-repr/structure/representation/ball-and-stick';
import { Download } from '../transforms/data';

export { DownloadAtomicStructure }
namespace DownloadAtomicStructure {
    export type Sources = 'pdbe-updated' | 'rcsb' | 'bcif-static' | 'url' | 'file'
    export type Source = ObtainStructureHelpers.MapParams<Sources, typeof ObtainStructureHelpers.ControlMap>
    export interface Params {
        source: Source
    }
}
namespace ObtainStructureHelpers {
    export const SourceOptions: [DownloadAtomicStructure.Sources, string][] = [
        ['pdbe-updated', 'PDBe Updated'],
        ['rcsb', 'RCSB'],
        ['bcif-static', 'BinaryCIF (static PDBe Updated)'],
        ['url', 'URL'],
        // TODO
        // ['file', 'File']
    ];
    export const ControlMap = {
        'pdbe-updated': PD.Text('1cbs', { label: 'Id' }),
        'rcsb': PD.Text('1tqn', { label: 'Id' }),
        'bcif-static': PD.Text('1tqn', { label: 'Id' }),
        'url': PD.Group({ url: PD.Text(''), isBinary: PD.Boolean(false) }, { isExpanded: true }),
        'file': PD.Group({ })
    }
    export function getControls(key: string) { return (ControlMap as any)[key]; }

    export type MapParams<P extends DownloadAtomicStructure.Sources, Map extends { [K in P]: PD.Any }> = P extends DownloadAtomicStructure.Sources ? PD.NamedParams<Map[P]['defaultValue'], P> : never

    export function getUrl(src: DownloadAtomicStructure.Source): Download.Params {
        switch (src.name) {
            case 'url': return src.params;
            case 'pdbe-updated': return { url: `https://www.ebi.ac.uk/pdbe/static/entry/${src.params.toLowerCase()}_updated.cif`, isBinary: false, label: `PDBe: ${src.params}` };
            case 'rcsb': return { url: `https://files.rcsb.org/download/${src.params.toUpperCase()}.cif`, isBinary: false, label: `RCSB: ${src.params}` };
            case 'bcif-static': return { url: `https://webchem.ncbr.muni.cz/ModelServer/static/bcif/${src.params.toLowerCase()}`, isBinary: true, label: `BinaryCIF: ${src.params}` };
            default: throw new Error(`${src.name} not supported.`);
        }
    }
}
const DownloadAtomicStructure = StateAction.create<PluginStateObject.Root, void, DownloadAtomicStructure.Params>({
    from: [PluginStateObject.Root],
    display: {
        name: 'Download Structure',
        description: 'Load a structure from PDBe and create its default Assembly and visual'
    },
    params: () => ({ source: PD.Mapped('bcif-static', ObtainStructureHelpers.SourceOptions, ObtainStructureHelpers.getControls) }),
    apply({ params, state }) {
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

        const url = ObtainStructureHelpers.getUrl(params.source);

        const root = b.toRoot()
            .apply(StateTransforms.Data.Download, url)
            .apply(StateTransforms.Data.ParseCif)
            .apply(StateTransforms.Model.TrajectoryFromMmCif, {})
            .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 })
            .apply(StateTransforms.Model.StructureAssemblyFromModel);

        root.apply(StateTransforms.Model.StructureComplexElement, { type: 'sequence' })
            .apply(StateTransforms.Representation.StructureRepresentation3D, { type: { name: 'cartoon', params: PD.getDefaultValues(CartoonParams) } });
        root.apply(StateTransforms.Model.StructureComplexElement, { type: 'ligands' })
            .apply(StateTransforms.Representation.StructureRepresentation3D, { type: { name: 'ball-and-stick', params: PD.getDefaultValues(BallAndStickParams) } });
        root.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' })
            .apply(StateTransforms.Representation.StructureRepresentation3D, { type: { name: 'ball-and-stick', params: { ...PD.getDefaultValues(BallAndStickParams), alpha: 0.51 } } });

        return state.update(root.getTree());
    }
});

export const CreateComplexRepresentation = StateAction.create<PluginStateObject.Molecule.Structure, void, {}>({
    from: [PluginStateObject.Molecule.Structure],
    display: {
        name: 'Create Complex',
        description: 'Split the structure into Sequence/Water/Ligands/... '
    },
    apply({ ref, state }) {
        const root = state.build().to(ref);
        root.apply(StateTransforms.Model.StructureComplexElement, { type: 'sequence' })
            .apply(StateTransforms.Representation.StructureRepresentation3D, { type: { name: 'cartoon', params: PD.getDefaultValues(CartoonParams) } });
        root.apply(StateTransforms.Model.StructureComplexElement, { type: 'ligands' })
            .apply(StateTransforms.Representation.StructureRepresentation3D, { type: { name: 'ball-and-stick', params: PD.getDefaultValues(BallAndStickParams) } });
        root.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' })
            .apply(StateTransforms.Representation.StructureRepresentation3D, { type: { name: 'ball-and-stick', params: { ...PD.getDefaultValues(BallAndStickParams), alpha: 0.51 } } });

        return state.update(root.getTree());
    }
});

export const UpdateTrajectory = StateAction.create<PluginStateObject.Root, void, { action: 'advance' | 'reset', by?: number }>({
    from: [],
    display: {
        name: 'Update Trajectory'
    },
    params: () => ({
        action: PD.Select('advance', [['advance', 'Advance'], ['reset', 'Reset']]),
        by: PD.Numeric(1, { min: -1, max: 1, step: 1 }, { isOptional: true })
    }),
    apply({ params, state }) {
        const models = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Model)
            .filter(c => c.transform.transformer === StateTransforms.Model.ModelFromTrajectory));

        const update = state.build();

        if (params.action === 'reset') {
            for (const m of models) {
                update.to(m.transform.ref).update(StateTransforms.Model.ModelFromTrajectory,
                    () => ({ modelIndex: 0}));
            }
        } else {
            for (const m of models) {
                const parent = StateSelection.findAncestorOfType(state.tree, state.cells, m.transform.ref, [PluginStateObject.Molecule.Trajectory]);
                if (!parent || !parent.obj) continue;
                const traj = parent.obj as PluginStateObject.Molecule.Trajectory;
                update.to(m.transform.ref).update(StateTransforms.Model.ModelFromTrajectory,
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