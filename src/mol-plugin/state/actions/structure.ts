/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from 'mol-plugin/context';
import { StateAction, StateBuilder, StateSelection, StateTransformer, State } from 'mol-state';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { Download } from '../transforms/data';
import { StructureRepresentation3DHelpers } from '../transforms/representation';
import { CustomModelProperties } from '../transforms/model';
import { DataFormatProvider } from './data-format';
import { FileInfo } from 'mol-util/file-info';
import { Task } from 'mol-task';

export const MmcifProvider: DataFormatProvider<any> = {
    label: 'mmCIF',
    description: 'mmCIF',
    stringExtensions: ['cif', 'mmcif', 'mcif'],
    binaryExtensions: ['bcif'],
    isApplicable: (info: FileInfo, data: Uint8Array) => {
        return info.ext === 'cif' || info.ext === 'mmcif' || info.ext === 'mcif' || info.ext === 'bcif'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, state: State) => {
        return Task.create('mmCIF default builder', async taskCtx => {
            const traj = createModelTree(data, 'cif');
            await state.updateTree(createStructureTree(ctx, traj, false)).runInContext(taskCtx)
        })
    }
}

export const PdbProvider: DataFormatProvider<any> = {
    label: 'PDB',
    description: 'PDB',
    stringExtensions: ['pdb', 'ent'],
    binaryExtensions: [],
    isApplicable: (info: FileInfo, data: Uint8Array) => {
        return info.ext === 'pdb' || info.ext === 'ent'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.String>, state: State) => {
        return Task.create('mmCIF default builder', async taskCtx => {
            const traj = createModelTree(data, 'pdb');
            await state.updateTree(createStructureTree(ctx, traj, false)).runInContext(taskCtx)
        })
    }
}

//

export { DownloadStructure };
type DownloadStructure = typeof DownloadStructure
const DownloadStructure = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Structure', description: 'Load a structure from the provided source and create its default Assembly and visual.' },
    params: {
        source: PD.MappedStatic('bcif-static', {
            'pdbe-updated': PD.Group({
                id: PD.Text('1cbs', { label: 'Id' }),
                supportProps: PD.Boolean(false)
            }, { isFlat: true }),
            'rcsb': PD.Group({
                id: PD.Text('1tqn', { label: 'Id' }),
                supportProps: PD.Boolean(false)
            }, { isFlat: true }),
            'bcif-static': PD.Group({
                id: PD.Text('1tqn', { label: 'Id' }),
                supportProps: PD.Boolean(false)
            }, { isFlat: true }),
            'url': PD.Group({
                url: PD.Text(''),
                format: PD.Select('cif', [['cif', 'CIF'], ['pdb', 'PDB']]),
                isBinary: PD.Boolean(false),
                supportProps: PD.Boolean(false)
            }, { isFlat: true })
        }, {
            options: [
                ['pdbe-updated', 'PDBe Updated'],
                ['rcsb', 'RCSB'],
                ['bcif-static', 'BinaryCIF (static PDBe Updated)'],
                ['url', 'URL']
            ]
        })
    }
})(({ params, state }, ctx: PluginContext) => {
    const b = state.build();
    const src = params.source;
    let downloadParams: StateTransformer.Params<Download>;

    switch (src.name) {
        case 'url':
            downloadParams = { url: src.params.url, isBinary: src.params.isBinary };
            break;
        case 'pdbe-updated':
            downloadParams = { url: `https://www.ebi.ac.uk/pdbe/static/entry/${src.params.id.toLowerCase()}_updated.cif`, isBinary: false, label: `PDBe: ${src.params.id}` };
            break;
        case 'rcsb':
            downloadParams = { url: `https://files.rcsb.org/download/${src.params.id.toUpperCase()}.cif`, isBinary: false, label: `RCSB: ${src.params.id}` };
            break;
        case 'bcif-static':
            downloadParams = { url: `https://webchem.ncbr.muni.cz/ModelServer/static/bcif/${src.params.id.toLowerCase()}`, isBinary: true, label: `BinaryCIF: ${src.params.id}` };
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    const data = b.toRoot().apply(StateTransforms.Data.Download, downloadParams, { props: { isGhost: true }});
    const traj = createModelTree(data, src.name === 'url' ? src.params.format : 'cif');
    return state.updateTree(createStructureTree(ctx, traj, params.source.params.supportProps));
});

function createModelTree(b: StateBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, format: 'pdb' | 'cif' = 'cif') {
    const parsed = format === 'cif'
        ? b.apply(StateTransforms.Data.ParseCif, void 0, { props: { isGhost: true }}).apply(StateTransforms.Model.TrajectoryFromMmCif, void 0, { props: { isGhost: true }})
        : b.apply(StateTransforms.Model.TrajectoryFromPDB, void 0, { props: { isGhost: true }});

    return parsed.apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 });
}

function createStructureTree(ctx: PluginContext, b: StateBuilder.To<PluginStateObject.Molecule.Model>, supportProps: boolean) {
    let root = b;
    if (supportProps) {
        root = root.apply(StateTransforms.Model.CustomModelProperties);
    }
    const structure = root.apply(StateTransforms.Model.StructureAssemblyFromModel);
    complexRepresentation(ctx, structure);

    return root;
}

function complexRepresentation(ctx: PluginContext, root: StateBuilder.To<PluginStateObject.Molecule.Structure>) {
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-sequence' })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'cartoon'));
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'atomic-het' })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick'));
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'water' })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'ball-and-stick', { alpha: 0.51 }));
    root.apply(StateTransforms.Model.StructureComplexElement, { type: 'spheres' })
        .apply(StateTransforms.Representation.StructureRepresentation3D,
            StructureRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'spacefill'));
}

export const CreateComplexRepresentation = StateAction.build({
    display: { name: 'Create Complex', description: 'Split the structure into Sequence/Water/Ligands/... ' },
    from: PluginStateObject.Molecule.Structure
})(({ ref, state }, ctx: PluginContext) => {
    const root = state.build().to(ref);
    complexRepresentation(ctx, root);
    return state.updateTree(root);
});

export const UpdateTrajectory = StateAction.build({
    display: { name: 'Update Trajectory' },
    params: {
        action: PD.Select<'advance' | 'reset'>('advance', [['advance', 'Advance'], ['reset', 'Reset']]),
        by: PD.makeOptional(PD.Numeric(1, { min: -1, max: 1, step: 1 }))
    }
})(({ params, state }) => {
    const models = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Model)
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

    return state.updateTree(update);
});

export const EnableModelCustomProps = StateAction.build({
    display: { name: 'Custom Properties', description: 'Enable the addition of custom properties to the model.' },
    from: PluginStateObject.Molecule.Model,
    params(a, ctx: PluginContext) {
        if (!a) return { properties: PD.MultiSelect([], [], { description: 'A list of property descriptor ids.' }) };
        return { properties: ctx.customModelProperties.getSelect(a.data) };
    },
    isApplicable(a, t, ctx: PluginContext) {
        return t.transformer !== CustomModelProperties;
    }
})(({ ref, params, state }, ctx: PluginContext) => {
    const root = state.build().to(ref).insert(CustomModelProperties, params);
    return state.updateTree(root);
});
