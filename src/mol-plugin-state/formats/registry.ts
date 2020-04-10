/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import msgpackDecode from '../../mol-io/common/msgpack/decode';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectRef } from '../../mol-state';
import { FileInfo } from '../../mol-util/file-info';
import { PluginStateObject } from '../objects';

export class DataFormatRegistry<Provider extends DataFormatProvider> {
    private _list: { name: string, provider: Provider }[] = []
    private _map = new Map<string, Provider>()
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

    constructor(buildInFormats: ReadonlyArray<readonly [string, Provider]>) {
        for (const [id, p] of buildInFormats) this.add(id, p);
    };

    private _clear() {
        this._extensions = undefined
        this._binaryExtensions = undefined
        this._options = undefined
    }

    add(name: string, provider: Provider) {
        this._clear()
        this._list.push({ name, provider })
        this._map.set(name, provider)
    }

    remove(name: string) {
        this._clear()
        this._list.splice(this._list.findIndex(e => e.name === name), 1)
        this._map.delete(name)
    }

    auto(info: FileInfo, dataStateObject: PluginStateObject.Data.Binary | PluginStateObject.Data.String) {
        for (let i = 0, il = this.list.length; i < il; ++i) {
            const { provider } = this._list[i];
            if (provider.isApplicable(info, dataStateObject.data)) return provider;
        }
        return;
    }

    get(name: string): Provider | undefined {
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

export interface DataFormatProvider<P = any, R = any> {
    label: string
    description: string
    stringExtensions: string[]
    binaryExtensions: string[]
    isApplicable(info: FileInfo, data: string | Uint8Array): boolean
    parse(plugin: PluginContext, data: StateObjectRef<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, params?: P): Promise<R>
}

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