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
import { CustomModelProperties, StructureSelection } from '../transforms/model';
import { DataFormatProvider } from './data-format';
import { FileInfo } from 'mol-util/file-info';
import { Task } from 'mol-task';
import { StructureElement } from 'mol-model/structure';

export const MmcifProvider: DataFormatProvider<any> = {
    label: 'mmCIF',
    description: 'mmCIF',
    stringExtensions: ['cif', 'mmcif', 'mcif'],
    binaryExtensions: ['bcif'],
    isApplicable: (info: FileInfo, data: Uint8Array | string) => {
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
    isApplicable: (info: FileInfo, data: string) => {
        return info.ext === 'pdb' || info.ext === 'ent'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.String>, state: State) => {
        return Task.create('PDB default builder', async taskCtx => {
            const traj = createModelTree(data, 'pdb');
            await state.updateTree(createStructureTree(ctx, traj, false)).runInContext(taskCtx)
        })
    }
}

export const GroProvider: DataFormatProvider<any> = {
    label: 'GRO',
    description: 'GRO',
    stringExtensions: ['gro'],
    binaryExtensions: [],
    isApplicable: (info: FileInfo, data: string) => {
        return info.ext === 'gro'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.String>, state: State) => {
        return Task.create('GRO default builder', async taskCtx => {
            const traj = createModelTree(data, 'gro');
            await state.updateTree(createStructureTree(ctx, traj, false)).runInContext(taskCtx)
        })
    }
}

//

const DownloadStructurePdbIdSourceOptions = PD.Group({
    supportProps: PD.asOptional(PD.Boolean(false)),
    asTrajectory: PD.asOptional(PD.Boolean(false, { description: 'Load all entries into a single trajectory.' }))
});

export { DownloadStructure };
type DownloadStructure = typeof DownloadStructure
const DownloadStructure = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Structure', description: 'Load a structure from the provided source and create its default Assembly and visual.' },
    params: {
        source: PD.MappedStatic('bcif-static', {
            'pdbe-updated': PD.Group({
                id: PD.Text('1cbs', { label: 'Id' }),
                options: DownloadStructurePdbIdSourceOptions
            }, { isFlat: true }),
            'rcsb': PD.Group({
                id: PD.Text('1tqn', { label: 'Id' }),
                options: DownloadStructurePdbIdSourceOptions
            }, { isFlat: true }),
            'bcif-static': PD.Group({
                id: PD.Text('1tqn', { label: 'Id' }),
                options: DownloadStructurePdbIdSourceOptions
            }, { isFlat: true }),
            'url': PD.Group({
                url: PD.Text(''),
                format: PD.Select('cif', [['cif', 'CIF'], ['pdb', 'PDB']]),
                isBinary: PD.Boolean(false),
                options: PD.Group({
                    supportProps: PD.asOptional(PD.Boolean(false))
                })
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
    let downloadParams: StateTransformer.Params<Download>[];
    let supportProps = false, asTrajectory = false;

    switch (src.name) {
        case 'url':
            downloadParams = [{ url: src.params.url, isBinary: src.params.isBinary }];
            supportProps = !!src.params.options.supportProps;
            break;
        case 'pdbe-updated':
            downloadParams = getDownloadParams(src.params.id, id => `https://www.ebi.ac.uk/pdbe/static/entry/${id.toLowerCase()}_updated.cif`, id => `PDBe: ${id}`, false);
            supportProps = !!src.params.options.supportProps;
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'rcsb':
            downloadParams = getDownloadParams(src.params.id, id => `https://files.rcsb.org/download/${id.toUpperCase()}.cif`, id => `RCSB: ${id}`, false);
            supportProps = !!src.params.options.supportProps;
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        case 'bcif-static':
            downloadParams = getDownloadParams(src.params.id, id => `https://webchem.ncbr.muni.cz/ModelServer/static/bcif/${id.toLowerCase()}`, id => `BinaryCIF: ${id}`, true);
            supportProps = !!src.params.options.supportProps;
            asTrajectory = !!src.params.options.asTrajectory;
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    if (downloadParams.length > 0 && asTrajectory) {
        const traj = createSingleTrajectoryModel(downloadParams, b);
        createStructureTree(ctx, traj, supportProps);
    } else {
        for (const download of downloadParams) {
            const data = b.toRoot().apply(StateTransforms.Data.Download, download, { props: { isGhost: true }});
            const traj = createModelTree(data, src.name === 'url' ? src.params.format : 'cif');
            createStructureTree(ctx, traj, supportProps)
        }
    }
    return state.updateTree(b);
});

function getDownloadParams(src: string, url: (id: string) => string, label: (id: string) => string, isBinary: boolean): StateTransformer.Params<Download>[] {
    const ids = src.split(',').map(id => id.trim()).filter(id => !!id && id.length >= 4);
    const ret: StateTransformer.Params<Download>[] = [];
    for (const id of ids) {
        ret.push({ url: url(id), isBinary, label: label(id) })
    }
    return ret;
}

function createSingleTrajectoryModel(sources: StateTransformer.Params<Download>[], b: StateBuilder.Root) {
    return b.toRoot()
        .apply(StateTransforms.Data.DownloadBlob, {
            sources: sources.map((src, i) => ({ id: '' + i, url: src.url, isBinary: src.isBinary })),
            maxConcurrency: 6
        }).apply(StateTransforms.Data.ParseBlob, {
            formats: sources.map((_, i) => ({ id: '' + i, format: 'cif' as 'cif' }))
        })
        .apply(StateTransforms.Model.TrajectoryFromBlob)
        .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 });
}

function createModelTree(b: StateBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, format: 'pdb' | 'cif' | 'gro' = 'cif') {
    let parsed: StateBuilder.To<PluginStateObject.Molecule.Trajectory>
    switch (format) {
        case 'cif':
            parsed = b.apply(StateTransforms.Data.ParseCif, void 0, { props: { isGhost: true }})
                .apply(StateTransforms.Model.TrajectoryFromMmCif, void 0, { props: { isGhost: true }})
            break
        case 'pdb':
            parsed = b.apply(StateTransforms.Model.TrajectoryFromPDB, void 0, { props: { isGhost: true }});
            break
        case 'gro':
            parsed = b.apply(StateTransforms.Model.TrajectoryFromGRO, void 0, { props: { isGhost: true }});
            break
        default:
            throw new Error('unsupported format')
    }

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
        by: PD.asOptional(PD.Numeric(1, { min: -1, max: 1, step: 1 }))
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

export const StructureFromSelection = StateAction.build({
    display: { name: 'Selection Structure', description: 'Create a new Structure from the current selection.' },
    from: PluginStateObject.Molecule.Structure,
    params: {
        label: PD.Text()
    }
    // isApplicable(a, t, ctx: PluginContext) {
    //     return t.transformer !== CustomModelProperties;
    // }
})(({ a, ref, params, state }, plugin: PluginContext) => {
    const sel = plugin.helpers.structureSelection.get(a.data);
    if (sel.kind === 'empty-loci') return Task.constant('', void 0);

    const query = StructureElement.Loci.toScriptExpression(sel);
    const root = state.build().to(ref).apply(StructureSelection, { query, label: params.label });
    return state.updateTree(root);
});


export const TestBlob = StateAction.build({
    display: { name: 'Test Blob'},
    from: PluginStateObject.Root
})(({ ref, state }, ctx: PluginContext) => {

    const ids = '5B6V,5B6W,5H2H,5H2I,5H2J,5B6X,5H2K,5H2L,5H2M,5B6Y,5H2N,5H2O,5H2P,5B6Z'.split(',').map(u => u.toLowerCase());

    const root = state.build().to(ref)
        .apply(StateTransforms.Data.DownloadBlob, {
            sources: ids.map(id => ({ id, url: `https://webchem.ncbr.muni.cz/ModelServer/static/bcif/${id}`, isBinary: true })),
            maxConcurrency: 4
        }).apply(StateTransforms.Data.ParseBlob, {
            formats: ids.map(id => ({ id, format: 'cif' as 'cif' }))
        })
        .apply(StateTransforms.Model.TrajectoryFromBlob)
        .apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 });
    createStructureTree(ctx, root, false);
    return state.updateTree(root);
});