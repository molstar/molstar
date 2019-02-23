/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeIsoValue } from 'mol-model/volume';
import { PluginContext } from 'mol-plugin/context';
import { State, StateAction, StateBuilder, StateObject, StateTransformer } from 'mol-state';
import { Task } from 'mol-task';
import { ColorNames } from 'mol-util/color/tables';
import { FileInfo, getFileInfo } from 'mol-util/file-info';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { Download } from '../transforms/data';
import { VolumeRepresentation3DHelpers } from '../transforms/representation';
import { VolumeStreaming } from 'mol-plugin/behavior/dynamic/volume';

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
    getDefaultBuilder(ctx: PluginContext, data: StateBuilder.To<D>, state?: State): Task<void>
}

const Ccp4Provider: DataFormatProvider<any> = {
    label: 'CCP4/MRC/BRIX',
    description: 'CCP4/MRC/BRIX',
    fileExtensions: ['ccp4', 'mrc', 'map'],
    isApplicable: (info: FileInfo, data: Uint8Array) => {
        return info.ext === 'ccp4' || info.ext === 'mrc' || info.ext === 'map'
    },
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.Binary>, state: State) => {
        return Task.create('CCP4/MRC/BRIX default builder', async taskCtx => {
            const tree = data.apply(StateTransforms.Data.ParseCcp4)
                .apply(StateTransforms.Volume.VolumeFromCcp4)
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
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.Binary>, state: State) => {
        return Task.create('DSN6/BRIX default builder', async taskCtx => {
            const tree = data.apply(StateTransforms.Data.ParseDsn6)
                .apply(StateTransforms.Volume.VolumeFromDsn6)
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
    getDefaultBuilder: (ctx: PluginContext, data: StateBuilder.To<PluginStateObject.Data.Binary>, state: State) => {
        return Task.create('DensityServer CIF default builder', async taskCtx => {
            const cifBuilder = data.apply(StateTransforms.Data.ParseCif)
            const cifStateObject = await state.updateTree(cifBuilder).runInContext(taskCtx)
            const b = state.build().to(cifBuilder.ref);
            const blocks = cifStateObject.data.blocks.slice(1); // zero block contains query meta-data
            let tree: StateBuilder.To<any>
            if (blocks.length === 1) {
                tree = b
                    .apply(StateTransforms.Volume.VolumeFromDensityServerCif, { blockHeader: blocks[0].header })
                    .apply(StateTransforms.Representation.VolumeRepresentation3D, VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface', { isoValue: VolumeIsoValue.relative(1.5), alpha: 0.3 }, 'uniform', { value: ColorNames.teal }))
            } else if (blocks.length === 2) {
                tree = b
                    .apply(StateTransforms.Volume.VolumeFromDensityServerCif, { blockHeader: blocks[0].header })
                    .apply(StateTransforms.Representation.VolumeRepresentation3D, VolumeRepresentation3DHelpers.getDefaultParamsStatic(ctx, 'isosurface', { isoValue: VolumeIsoValue.relative(1.5), alpha: 0.3 }, 'uniform', { value: ColorNames.blue }))
                const vol = tree.to(cifBuilder.ref)
                    .apply(StateTransforms.Volume.VolumeFromDensityServerCif, { blockHeader: blocks[1].header })
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

    // Alternative for more complex states where the builder is not a simple StateBuilder.To<>:
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
                'pdbe-emd-ds': PD.Group({
                    id: PD.Text('emd-8004', { label: 'Id' }),
                    detail: PD.Numeric(3, { min: 0, max: 10, step: 1 }, { label: 'Detail' }),
                }, { isFlat: true }),
                'pdbe-xray-ds': PD.Group({
                    id: PD.Text('1tqn', { label: 'Id' }),
                    detail: PD.Numeric(3, { min: 0, max: 10, step: 1 }, { label: 'Detail' }),
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
                    ['pdbe-emd-ds', 'PDBe EMD Density Server'],
                    ['pdbe-xray-ds', 'PDBe X-ray Density Server'],
                    ['rcsb', 'RCSB X-ray maps'],
                    ['url', 'URL']
                ]
            })
        }
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Download Density', async taskCtx => {
    const src = params.source;
    let downloadParams: StateTransformer.Params<Download>;
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
        case 'pdbe-emd-ds':
            downloadParams = {
                url: `https://www.ebi.ac.uk/pdbe/densities/emd/${src.params.id.toLowerCase()}/cell?detail=${src.params.detail}`,
                isBinary: true,
                label: `PDBe EMD Density Server: ${src.params.id}`
            };
            break;
        case 'pdbe-xray-ds':
            downloadParams = {
                url: `https://www.ebi.ac.uk/pdbe/densities/x-ray/${src.params.id.toLowerCase()}/cell?detail=${src.params.detail}`,
                isBinary: true,
                label: `PDBe X-ray Density Server: ${src.params.id}`
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
        case 'pdbe-emd-ds':
        case 'pdbe-xray-ds':
            provider = ctx.dataFormat.registry.get('dscif')
            break;
        case 'rcsb':
            provider = ctx.dataFormat.registry.get('dsn6')
            break;
        default: throw new Error(`${(src as any).name} not supported.`);
    }

    const b = state.build().to(data.ref);
    await provider.getDefaultBuilder(ctx, b, state).runInContext(taskCtx)
}));

export const InitVolumeStreaming = StateAction.build({
    display: { name: 'Volume Streaming' },
    from: PluginStateObject.Molecule.Model,
    params: VolumeStreaming.Params
})(({ ref, state, params }, ctx: PluginContext) => {
    // TODO: specify simpler params
    // TODO: try to determine if the input is x-ray or emd (in params provider)
    // TODO: for EMD, use PDBe API to determine controur level https://github.com/dsehnal/LiteMol/blob/master/src/Viewer/Extensions/DensityStreaming/Entity.ts#L168
    // TODO: custom react view for this and the VolumeStreamingBehavior transformer

    const root = state.build().to(ref)
        .apply(StateTransforms.Volume.VolumeStreamingBehavior, params);

    root.apply(StateTransforms.Volume.VolumeStreamingVisual, { channel: '2FO-FC', level: '2fo-fc' });
    root.apply(StateTransforms.Volume.VolumeStreamingVisual, { channel: 'FO-FC', level: 'fo-fc(+ve)' });
    root.apply(StateTransforms.Volume.VolumeStreamingVisual, { channel: 'FO-FC', level: 'fo-fc(-ve)' });

    return state.updateTree(root);
});