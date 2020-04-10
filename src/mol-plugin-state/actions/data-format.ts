/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import msgpackDecode from '../../mol-io/common/msgpack/decode';
import { PluginContext } from '../../mol-plugin/context';
import { State, StateAction, StateObjectRef } from '../../mol-state';
import { Task } from '../../mol-task';
import { FileInfo, getFileInfo } from '../../mol-util/file-info';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateObject } from '../objects';
import { PlyProvider } from './shape';
import { DcdProvider, GroProvider, MmcifProvider, PdbProvider, Provider3dg, PsfProvider, MolProvider, CifCoreProvider } from './structure';
import { Ccp4Provider, DscifProvider, Dsn6Provider } from './volume';

export class DataFormatRegistry<D extends PluginStateObject.Data.Binary | PluginStateObject.Data.String> {
    private _list: { name: string, provider: DataFormatProvider<D> }[] = []
    private _map = new Map<string, DataFormatProvider<D>>()
    private _extensions: Set<string> | undefined = undefined
    private _binaryExtensions: Set<string> | undefined = undefined
    private _options: [string, string][] | undefined = undefined

    get types(): [string, string][] {
        return this._list.map(e => [e.name, e.provider.label] as [string, string]);
    }

    get extensions() {
        if (this._extensions) return this._extensions
        const extensions = new Set<string>()
        this._list.forEach(({ provider }) => {
            provider.stringExtensions.forEach(ext => extensions.add(ext))
            provider.binaryExtensions.forEach(ext => extensions.add(ext))
        })
        this._extensions = extensions
        return extensions
    }

    get binaryExtensions() {
        if (this._binaryExtensions) return this._binaryExtensions
        const binaryExtensions = new Set<string>()
        this._list.forEach(({ provider }) => provider.binaryExtensions.forEach(ext => binaryExtensions.add(ext)))
        this._binaryExtensions = binaryExtensions
        return binaryExtensions
    }

    get options() {
        if (this._options) return this._options
        const options: [string, string][] = [['auto', 'Automatic']]
        this._list.forEach(({ name, provider }) => options.push([ name, provider.label ]))
        this._options = options
        return options
    }

    constructor() {
        this.add('3dg', Provider3dg)
        this.add('ccp4', Ccp4Provider)
        this.add('cifCore', CifCoreProvider)
        this.add('dcd', DcdProvider)
        this.add('dscif', DscifProvider)
        this.add('dsn6', Dsn6Provider)
        this.add('gro', GroProvider)
        this.add('mol', MolProvider)
        this.add('mmcif', MmcifProvider)
        this.add('pdb', PdbProvider)
        this.add('ply', PlyProvider)
        this.add('psf', PsfProvider)
    };

    private _clear() {
        this._extensions = undefined
        this._binaryExtensions = undefined
        this._options = undefined
    }

    add(name: string, provider: DataFormatProvider<D>) {
        this._clear()
        this._list.push({ name, provider })
        this._map.set(name, provider)
    }

    remove(name: string) {
        this._clear()
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

export type DataFormatBuilderOptions = { visuals: boolean }

export interface DataFormatProvider<D extends PluginStateObject.Data.Binary | PluginStateObject.Data.String> {
    label: string
    description: string
    stringExtensions: string[]
    binaryExtensions: string[]
    isApplicable(info: FileInfo, data: string | Uint8Array): boolean
    getDefaultBuilder(ctx: PluginContext, data: StateObjectRef<D>, options: DataFormatBuilderOptions, state: State): Task<void>
}

//

export const OpenFiles = StateAction.build({
    display: { name: 'Open Files', description: 'Load one or more files and optionally create default visuals' },
    from: PluginStateObject.Root,
    params: (a, ctx: PluginContext) => {
        const { extensions, options } = ctx.dataFormat.registry
        return {
            files: PD.FileList({ accept: Array.from(extensions.values()).map(e => `.${e}`).join(',') + ',.gz,.zip', multiple: true }),
            format: PD.Select('auto', options),
            visuals: PD.Boolean(true, { description: 'Add default visuals' }),
        }
    }
})(({ params, state }, plugin: PluginContext) => Task.create('Open Files', async taskCtx => {
    await state.transaction(async () => {
        if (params.files === null) {
            plugin.log.error('No file(s) selected')
            return
        }
        for (let i = 0, il = params.files.length; i < il; ++i) {
            try {
                const file = params.files[i]
                const info = getFileInfo(file)
                const isBinary = plugin.dataFormat.registry.binaryExtensions.has(info.ext)
                const { data } = await plugin.builders.data.readFile({ file, isBinary });
                const provider = params.format === 'auto'
                    ? plugin.dataFormat.registry.auto(info, data.cell?.obj!)
                    : plugin.dataFormat.registry.get(params.format)

                // need to await so that the enclosing Task finishes after the update is done.
                await provider.getDefaultBuilder(plugin, data, { visuals: params.visuals }, state).runInContext(taskCtx)
            } catch (e) {
                plugin.log.error(e)
            }
        }
    }).runInContext(taskCtx);
}));

//

type cifVariants = 'dscif' | 'coreCif' | -1
export function guessCifVariant(info: FileInfo, data: Uint8Array | string): cifVariants {
    if (info.ext === 'bcif') {
        try {
            // TODO: find a way to run msgpackDecode only once
            //      now it is run twice, here and during file parsing
            if (msgpackDecode(data as Uint8Array).encoder.startsWith('VolumeServer')) return 'dscif'
        } catch { }
    } else if (info.ext === 'cif') {
        const str = data as string
        if (str.startsWith('data_SERVER\n#\n_density_server_result')) return 'dscif'
        if (str.includes('atom_site_fract_x') || str.includes('atom_site.fract_x')) return 'coreCif'
    }
    return -1
}