/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from 'mol-plugin/context';
import { StateTree, Transformer, StateObject, State } from 'mol-state';
import { StateAction } from 'mol-state/action';
import { StateSelection } from 'mol-state/state/selection';
import { StateTreeBuilder } from 'mol-state/tree/builder';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { Download } from '../transforms/data';
import { StructureRepresentation3DHelpers, VolumeRepresentation3DHelpers } from '../transforms/representation';
import { getFileInfo, FileInfo } from 'mol-util/file-info';
import { Task } from 'mol-task';
import { ColorNames } from 'mol-util/color/tables';
import { VolumeIsoValue } from 'mol-model/volume';

// TODO: "structure/volume parser provider"

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
    let downloadParams: Transformer.Params<Download>;

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

    const data = b.toRoot().apply(StateTransforms.Data.Download, downloadParams);
    const traj = createModelTree(data, src.name === 'url' ? src.params.format : 'cif');
    return state.updateTree(createStructureTree(ctx, traj, params.source.params.supportProps));
});

export const OpenStructure = StateAction.build({
    display: { name: 'Open Structure', description: 'Load a structure from file and create its default Assembly and visual' },
    from: PluginStateObject.Root,
    params: { file: PD.File({ accept: '.cif,.bcif' }) }
})(({ params, state }, ctx: PluginContext) => {
    const b = state.build();
    const data = b.toRoot().apply(StateTransforms.Data.ReadFile, { file: params.file, isBinary: /\.bcif$/i.test(params.file.name) });
    const traj = createModelTree(data, 'cif');
    return state.updateTree(createStructureTree(ctx, traj, false));
});

function createModelTree(b: StateTreeBuilder.To<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, format: 'pdb' | 'cif' = 'cif') {
    const parsed = format === 'cif'
        ? b.apply(StateTransforms.Data.ParseCif).apply(StateTransforms.Model.TrajectoryFromMmCif)
        : b.apply(StateTransforms.Model.TrajectoryFromPDB);

    return parsed.apply(StateTransforms.Model.ModelFromTrajectory, { modelIndex: 0 });
}

function createStructureTree(ctx: PluginContext, b: StateTreeBuilder.To<PluginStateObject.Molecule.Model>, supportProps: boolean): StateTree {
    let root = b;
    if (supportProps) {
        root = root.apply(StateTransforms.Model.CustomModelProperties);
    }
    const structure = root.apply(StateTransforms.Model.StructureAssemblyFromModel);
    complexRepresentation(ctx, structure);

    return root.getTree();
}

function complexRepresentation(ctx: PluginContext, root: StateTreeBuilder.To<PluginStateObject.Molecule.Structure>) {
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
    return state.updateTree(root.getTree());
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

//

export class DataFormatRegistry<D extends PluginStateObject.Data.Binary | PluginStateObject.Data.String, M extends StateObject> {
    private _list: { name: string, provider: DataFormatProvider<D> }[] = []
    private _map = new Map<string, DataFormatProvider<D>>()

    get default() { return this._list[0]; }
    get types(): [string, string][] {
        return this._list.map(e => [e.name, e.provider.label] as [string, string]);
    }

    constructor() {
        this.add('ccp4', Ccp4Provider)
        this.add('dsn6', Dsn6Provider)
        this.add('dscif', DscifProvider)
    };

    add(name: string, provider: DataFormatProvider<D>) {
        this._list.push({ name, provider })
        this._map.set(name, provider)
    }

    remove(name: string) {
        this._list.splice(this._list.findIndex(e => e.name === name), 1)
        this._map.delete(name)
    }

    auto(info: FileInfo, dataStateObject: D) {
        for (let i = 0, il = this.list.length; i < il; ++i) {
            const { provider } = this._list[i]
            if (provider.isApplicable(info, dataStateObject.data)) return provider
        }
        throw new Error('no compatible data format provider available')
    }

    get(name: string): DataFormatProvider<D> {
        if (this._map.has(name)) {
            return this._map.get(name)!
        } else {
            throw new Error(`unknown data format name '${name}'`)
        }
    }

    get list() {
        return this._list
    }
}

interface DataFormatProvider<D extends PluginStateObject.Data.Binary | PluginStateObject.Data.String> {
    label: string
    description: string
    fileExtensions: string[]
    isApplicable(info: FileInfo, data: string | Uint8Array): boolean
    getDefaultBuilder(ctx: PluginContext, data: StateTreeBuilder.To<D>, state?: State): Task<void>
}

const Ccp4Provider: DataFormatProvider<any> = {
    label: 'CCP4/MRC/BRIX',
    description: 'CCP4/MRC/BRIX',
    fileExtensions: ['ccp4', 'mrc', 'map'],
    isApplicable: (info: FileInfo, data: Uint8Array) => {
        return info.ext === 'ccp4' || info.ext === 'mrc' || info.ext === 'map'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateTreeBuilder.To<PluginStateObject.Data.Binary>, state: State) => {
        return Task.create('CCP4/MRC/BRIX default builder', async taskCtx => {
            const tree = data.apply(StateTransforms.Data.ParseCcp4)
                .apply(StateTransforms.Model.VolumeFromCcp4)
                .apply(StateTransforms.Representation.VolumeRepresentation3D)
            await state.updateTree(tree).runInContext(taskCtx)
        })
    }
}

const Dsn6Provider: DataFormatProvider<any> = {
    label: 'DSN6/BRIX',
    description: 'DSN6/BRIX',
    fileExtensions: ['dsn6', 'brix'],
    isApplicable: (info: FileInfo, data: Uint8Array) => {
        return info.ext === 'dsn6' || info.ext === 'brix'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateTreeBuilder.To<PluginStateObject.Data.Binary>, state: State) => {
        return Task.create('DSN6/BRIX default builder', async taskCtx => {
            const tree = data.apply(StateTransforms.Data.ParseDsn6)
                .apply(StateTransforms.Model.VolumeFromDsn6)
                .apply(StateTransforms.Representation.VolumeRepresentation3D)
            await state.updateTree(tree).runInContext(taskCtx)
        })
    }
}

const DscifProvider: DataFormatProvider<any> = {
    label: 'DensityServer CIF',
    description: 'DensityServer CIF',
    fileExtensions: ['cif'],
    isApplicable: (info: FileInfo, data: Uint8Array) => {
        return info.ext === 'cif'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateTreeBuilder.To<PluginStateObject.Data.Binary>, state: State) => {
        return Task.create('DensityServer CIF default builder', async taskCtx => {
            const cifBuilder = data.apply(StateTransforms.Data.ParseCif)
            const cifStateObject = await state.updateTree(cifBuilder).runInContext(taskCtx)
            const b = state.build().to(cifBuilder.ref);
            const blocks = cifStateObject.data.blocks.slice(1); // zero block contains query meta-data
            let tree: StateTreeBuilder.To<any>
            if (blocks.length === 1) {
                tree = b
                    .apply(StateTransforms.Model.VolumeFromDensityServerCif, { blockHeader: blocks[0].header })
                    .apply(StateTransforms.Representation.VolumeRepresentation3D, VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface', { isoValue: VolumeIsoValue.relative(1.5), alpha: 0.3 }))
            } else if (blocks.length === 2) {
                tree = b
                    .apply(StateTransforms.Model.VolumeFromDensityServerCif, { blockHeader: blocks[0].header })
                    .apply(StateTransforms.Representation.VolumeRepresentation3D, VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface', { isoValue: VolumeIsoValue.relative(1.5), alpha: 0.3 }, 'uniform', { value: ColorNames.blue }))
                const vol = tree.to(cifBuilder.ref)
                    .apply(StateTransforms.Model.VolumeFromDensityServerCif, { blockHeader: blocks[1].header })
                const posParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface', { isoValue: VolumeIsoValue.relative(3), alpha: 0.3 }, 'uniform', { value: ColorNames.green })
                tree = vol.apply(StateTransforms.Representation.VolumeRepresentation3D, posParams)
                const negParams = VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface', { isoValue: VolumeIsoValue.relative(-3), alpha: 0.3 }, 'uniform', { value: ColorNames.red })
                tree = tree.to(vol.ref).apply(StateTransforms.Representation.VolumeRepresentation3D, negParams)
            } else {
                throw new Error('unknown number of blocks')
            }

            await state.updateTree(tree).runInContext(taskCtx);
        })
    }
}

//

function getDataFormatExtensionsOptions(dataFormatRegistry: DataFormatRegistry<any, any>) {
    const extensions: string[] = []
    const options: [string, string][] = [['auto', 'Automatic']]
    dataFormatRegistry.list.forEach(({ name, provider }) => {
        extensions.push(...provider.fileExtensions)
        options.push([ name, provider.label ])
    })
    return { extensions, options }
}

export const OpenVolume = StateAction.build({
    display: { name: 'Open Volume', description: 'Load a volume from file and create its default visual' },
    from: PluginStateObject.Root,
    params: (a, ctx: PluginContext) => {
        const { extensions, options } = getDataFormatExtensionsOptions(ctx.dataFormat.registry)
        return {
            file: PD.File({ accept: extensions.map(e => `.${e}`).join(',')}),
            format: PD.Select('auto', options),
            isBinary: PD.Boolean(true), // TOOD should take selected format into account
        }
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Open Volume', async taskCtx => {
    const data = state.build().toRoot().apply(StateTransforms.Data.ReadFile, { file: params.file, isBinary: params.isBinary });
    const dataStateObject = await state.updateTree(data).runInContext(taskCtx);

    // Alternative for more complex states where the builder is not a simple StateTreeBuilder.To<>:
    /*
    const dataRef = dataTree.ref;
    await state.updateTree(dataTree).runInContext(taskCtx);
    const dataCell = state.select(dataRef)[0];
    */

    const provider = params.format === 'auto' ? ctx.dataFormat.registry.auto(getFileInfo(params.file), dataStateObject) : ctx.dataFormat.registry.get(params.format)
    const b = state.build().to(data.ref);
    // need to await the 2nd update the so that the enclosing Task finishes after the update is done.
    await provider.getDefaultBuilder(ctx, b, state).runInContext(taskCtx)
}));

export { DownloadDensity };
type DownloadDensity = typeof DownloadDensity
const DownloadDensity = StateAction.build({
    from: PluginStateObject.Root,
    display: { name: 'Download Density', description: 'Load a density from the provided source and create its default visual.' },
    params: (a, ctx: PluginContext) => {
        const { options } = getDataFormatExtensionsOptions(ctx.dataFormat.registry)
        return {
            source: PD.MappedStatic('rcsb', {
                'pdbe': PD.Group({
                    id: PD.Text('1tqn', { label: 'Id' }),
                    type: PD.Select('2fofc', [['2fofc', '2Fo-Fc'], ['fofc', 'Fo-Fc']]),
                }, { isFlat: true }),
                'rcsb': PD.Group({
                    id: PD.Text('1tqn', { label: 'Id' }),
                    type: PD.Select('2fofc', [['2fofc', '2Fo-Fc'], ['fofc', 'Fo-Fc']]),
                }, { isFlat: true }),
                'url': PD.Group({
                    url: PD.Text(''),
                    isBinary: PD.Boolean(false),
                    format: PD.Select('auto', options),
                }, { isFlat: true })
            }, {
                options: [
                    ['pdbe', 'PDBe X-ray maps'],
                    ['rcsb', 'RCSB X-ray maps'],
                    ['url', 'URL']
                ]
            })
        }
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Download Density', async taskCtx => {
    const src = params.source;
    let downloadParams: Transformer.Params<Download>;
    let provider: DataFormatProvider<any>

    switch (src.name) {
        case 'url':
            downloadParams = src.params;
            break;
        case 'pdbe':
            downloadParams = {
                url: src.params.type === '2fofc'
                    ? `http://www.ebi.ac.uk/pdbe/coordinates/files/${src.params.id.toLowerCase()}.ccp4`
                    : `http://www.ebi.ac.uk/pdbe/coordinates/files/${src.params.id.toLowerCase()}_diff.ccp4`,
                isBinary: true,
                label: `PDBe X-ray map: ${src.params.id}`
            };
            break;
        case 'rcsb':
            downloadParams = {
                url: src.params.type === '2fofc'
                    ? `https://edmaps.rcsb.org/maps/${src.params.id.toLowerCase()}_2fofc.dsn6`
                    : `https://edmaps.rcsb.org/maps/${src.params.id.toLowerCase()}_fofc.dsn6`,
                isBinary: true,
                label: `RCSB X-ray map: ${src.params.id}`
            };
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    const data = state.build().toRoot().apply(StateTransforms.Data.Download, downloadParams);
    const dataStateObject = await state.updateTree(data).runInContext(taskCtx);

    switch (src.name) {
        case 'url':
            downloadParams = src.params;
            provider = src.params.format === 'auto' ? ctx.dataFormat.registry.auto(getFileInfo(downloadParams.url), dataStateObject) : ctx.dataFormat.registry.get(src.params.format)
            break;
        case 'pdbe':
            provider = ctx.dataFormat.registry.get('ccp4')
            break;
        case 'rcsb':
            provider = ctx.dataFormat.registry.get('dsn6')
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    const b = state.build().to(data.ref);
    await provider.getDefaultBuilder(ctx, b, state).runInContext(taskCtx)
}));