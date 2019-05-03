/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from 'mol-plugin/context';
import { State, StateBuilder, StateAction } from 'mol-state';
import { Task } from 'mol-task';
import { FileInfo, getFileInfo } from 'mol-util/file-info';
import { PluginStateObject } from '../objects';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Ccp4Provider, Dsn6Provider, DscifProvider } from './volume';
import { StateTransforms } from '../transforms';
import { MmcifProvider, PdbProvider, GroProvider } from './structure';
import msgpackDecode from 'mol-io/common/msgpack/decode'
import { PlyProvider } from './shape';

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
        this.add('ccp4', Ccp4Provider)
        this.add('dscif', DscifProvider)
        this.add('dsn6', Dsn6Provider)
        this.add('gro', GroProvider)
        this.add('mmcif', MmcifProvider)
        this.add('pdb', PdbProvider)
        this.add('ply', PlyProvider)
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

export interface DataFormatProvider<D extends PluginStateObject.Data.Binary | PluginStateObject.Data.String> {
    label: string
    description: string
    stringExtensions: string[]
    binaryExtensions: string[]
    isApplicable(info: FileInfo, data: string | Uint8Array): boolean
    getDefaultBuilder(ctx: PluginContext, data: StateBuilder.To<D>, state?: State): Task<void>
}

//

export const OpenFile = StateAction.build({
    display: { name: 'Open File', description: 'Load a file and create its default visual' },
    from: PluginStateObject.Root,
    params: (a, ctx: PluginContext) => {
        const { extensions, options } = ctx.dataFormat.registry
        return {
            file: PD.File({ accept: Array.from(extensions).map(e => `.${e}`).join(',')}),
            format: PD.Select('auto', options),
        }
    }
})(({ params, state }, ctx: PluginContext) => Task.create('Open File', async taskCtx => {
    const info = getFileInfo(params.file)
    const data = state.build().toRoot().apply(StateTransforms.Data.ReadFile, { file: params.file, isBinary: ctx.dataFormat.registry.binaryExtensions.has(info.ext) });
    const dataStateObject = await state.updateTree(data).runInContext(taskCtx);

    // Alternative for more complex states where the builder is not a simple StateBuilder.To<>:
    /*
    const dataRef = dataTree.ref;
    await state.updateTree(dataTree).runInContext(taskCtx);
    const dataCell = state.select(dataRef)[0];
    */

    // const data = b.toRoot().apply(StateTransforms.Data.ReadFile, { file: params.file, isBinary: /\.bcif$/i.test(params.file.name) });

    const provider = params.format === 'auto' ? ctx.dataFormat.registry.auto(info, dataStateObject) : ctx.dataFormat.registry.get(params.format)
    const b = state.build().to(data.ref);
    // need to await the 2nd update the so that the enclosing Task finishes after the update is done.
    await provider.getDefaultBuilder(ctx, b, state).runInContext(taskCtx)
}));

//

type cifVariants = 'dscif' | -1
export function guessCifVariant(info: FileInfo, data: Uint8Array | string): cifVariants {
    if (info.ext === 'bcif') {
        try {
            // TODO find a way to run msgpackDecode only once
            //      now it is run twice, here and during file parsing
            if (msgpackDecode(data as Uint8Array).encoder.startsWith('VolumeServer')) return 'dscif'
        } catch { }
    } else if (info.ext === 'cif') {
        if ((data as string).startsWith('data_SERVER\n#\n_density_server_result')) return 'dscif'
    }
    return -1
}