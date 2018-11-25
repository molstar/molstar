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
import { StateTree, Transformer } from 'mol-state';
import { StateTreeBuilder } from 'mol-state/tree/builder';
import { PolymerIdColorThemeParams } from 'mol-theme/color/polymer-id';
import { UniformSizeThemeParams } from 'mol-theme/size/uniform';
import { ElementSymbolColorThemeParams } from 'mol-theme/color/element-symbol';

// TODO: "structure parser provider"

export { DownloadStructure }
type DownloadStructure = typeof DownloadStructure
const DownloadStructure = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Structure', description: 'Load a structure from the provided source and create its default Assembly and visual.' },
    params: {
        source: PD.MappedStatic('bcif-static', {
            'pdbe-updated': PD.Text('1cbs', { label: 'Id' }),
            'rcsb': PD.Text('1tqn', { label: 'Id' }),
            'bcif-static': PD.Text('1tqn', { label: 'Id' }),
            'url': PD.Group({ url: PD.Text(''), isBinary: PD.Boolean(false) }, { isExpanded: true })
        }, {
            options: [
                ['pdbe-updated', 'PDBe Updated'],
                ['rcsb', 'RCSB'],
                ['bcif-static', 'BinaryCIF (static PDBe Updated)'],
                ['url', 'URL']
            ]
        })
    }
})(({ params, state }) => {
    const b = state.build();
    const src = params.source;
    let url: Transformer.Params<Download>;

    switch (src.name) {
        case 'url':
            url = src.params;
            break;
        case 'pdbe-updated':
            url = { url: `https://www.ebi.ac.uk/pdbe/static/entry/${src.params.toLowerCase()}_updated.cif`, isBinary: false, label: `PDBe: ${src.params}` };
            break;
        case 'rcsb':
            url = { url: `https://files.rcsb.org/download/${src.params.toUpperCase()}.cif`, isBinary: false, label: `RCSB: ${src.params}` };
            break;
        case 'bcif-static':
            url = { url: `https://webchem.ncbr.muni.cz/ModelServer/static/bcif/${src.params.toLowerCase()}`, isBinary: true, label: `BinaryCIF: ${src.params}` };
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    const data = b.toRoot().apply(StateTransforms.Data.Download, url);
    return state.update(createStructureTree(data));
});

export const OpenStructure = StateAction.build({
    display: { name: 'Open Structure', description: 'Load a structure from file and create its default Assembly and visual' },
    from: PluginStateObject.Root,
    params: { file: PD.File({ accept: '.cif,.bcif' }) }
})(({ params, state }) => {
    const b = state.build();
    const data = b.toRoot().apply(StateTransforms.Data.ReadFile, { file: params.file, isBinary: /\.bcif$/i.test(params.file.name) });
    return state.update(createStructureTree(data));
});

function createStructureTree(b: StateTreeBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>): StateTree {
    const root = b
        .apply(StateTransforms.Data.ParseCif)
        .apply(StateTransforms.Model.TrajectoryFromMmCif, {})
        .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 })
        .apply(StateTransforms.Model.StructureAssemblyFromModel);

    complexRepresentation(root);

    return root.getTree();
}

function complexRepresentation(root: StateTreeBuilder.To<PluginStateObject.Molecule.Structure>) {
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-sequence' })
        .apply(StateTransforms.Representation.StructureRepresentation3D, {
            type: { name: 'cartoon', params: PD.getDefaultValues(CartoonParams) },
            colorTheme: { name: 'polymer-id', params: PD.getDefaultValues(PolymerIdColorThemeParams) },
            sizeTheme: { name: 'uniform', params: PD.getDefaultValues(UniformSizeThemeParams) },
        });
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-het' })
        .apply(StateTransforms.Representation.StructureRepresentation3D, {
            type: { name: 'ball-and-stick', params: PD.getDefaultValues(BallAndStickParams) },
            colorTheme: { name: 'element-symbol', params: PD.getDefaultValues(ElementSymbolColorThemeParams) },
            sizeTheme: { name: 'uniform', params: PD.getDefaultValues(UniformSizeThemeParams) },
        });
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' })
        .apply(StateTransforms.Representation.StructureRepresentation3D, {
            type: { name: 'ball-and-stick', params: { ...PD.getDefaultValues(BallAndStickParams), alpha: 0.51 } },
            colorTheme: { name: 'element-symbol', params: PD.getDefaultValues(ElementSymbolColorThemeParams) },
            sizeTheme: { name: 'uniform', params: PD.getDefaultValues(UniformSizeThemeParams) },
        })
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'spheres' });
    // TODO: create spheres visual
}

export const CreateComplexRepresentation = StateAction.build({
    display: { name: 'Create Complex', description: 'Split the structure into Sequence/Water/Ligands/... ' },
    from: PluginStateObject.Molecule.Structure
})(({ ref, state }) => {
    const root = state.build().to(ref);
    complexRepresentation(root);
    return state.update(root.getTree());
});

export const UpdateTrajectory = StateAction.build({
    display: { name: 'Update Trajectory' },
    params: {
        action: PD.Select<'advance' | 'reset'>('advance', [['advance', 'Advance'], ['reset', 'Reset']]),
        by: PD.makeOptional(PD.Numeric(1, { min: -1, max: 1, step: 1 }))
    }
})(({ params, state }) => {
    const models = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Model)
        .filter(c => c.transform.transformer === StateTransforms.Model.ModelFromTrajectory));

    const update = state.build();

    if (params.action === 'reset') {
        for (const m of models) {
            update.to(m.transform.ref).update(StateTransforms.Model.ModelFromTrajectory,
                () => ({ modelIndex: 0 }));
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
});