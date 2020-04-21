/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FileInfo } from '../../mol-util/file-info';
import { PluginStateObject } from '../objects';
import { DataFormatProvider } from './provider';
import { BuiltInTrajectoryFormats } from './trajectory';
import { BuiltInVolumeFormats } from './volume';
import { BuiltInShapeFormats } from './shape';
import { BuiltInStructureFormats } from './structure';

export class DataFormatRegistry {
    private _list: { name: string, provider: DataFormatProvider }[] = []
    private _map = new Map<string, DataFormatProvider>()
    private _extensions: Set<string> | undefined = undefined
    private _binaryExtensions: Set<string> | undefined = undefined
    private _options: [string, string, string][] | undefined = undefined

    get types(): [string, string][] {
        return this._list.map(e => [e.name, e.provider.label] as [string, string]);
    }

    get extensions() {
        if (this._extensions) return this._extensions;
        const extensions = new Set<string>();
        this._list.forEach(({ provider }) => {
            provider.stringExtensions?.forEach(ext => extensions.add(ext));
            provider.binaryExtensions?.forEach(ext => extensions.add(ext));
        });
        this._extensions = extensions;
        return extensions;
    }

    get binaryExtensions() {
        if (this._binaryExtensions) return this._binaryExtensions;
        const binaryExtensions = new Set<string>();
        this._list.forEach(({ provider }) => provider.binaryExtensions?.forEach(ext => binaryExtensions.add(ext)));
        this._binaryExtensions = binaryExtensions;
        return binaryExtensions;
    }

    get options() {
        if (this._options) return this._options;
        const options: [string, string, string][] = [];
        this._list.forEach(({ name, provider }) => options.push([ name, provider.label, provider.category || '' ]));
        this._options = options;
        return options;
    }

    constructor() {
        for (const [id, p] of BuiltInVolumeFormats) this.add(id, p);
        for (const [id, p] of BuiltInStructureFormats) this.add(id, p);
        for (const [id, p] of BuiltInShapeFormats) this.add(id, p);
        for (const [id, p] of BuiltInTrajectoryFormats) this.add(id, p);
    };

    private _clear() {
        this._extensions = undefined;
        this._binaryExtensions = undefined;
        this._options = undefined;
    }

    add(name: string, provider: DataFormatProvider) {
        this._clear();
        this._list.push({ name, provider });
        this._map.set(name, provider);
    }

    remove(name: string) {
        this._clear();
        this._list.splice(this._list.findIndex(e => e.name === name), 1);
        this._map.delete(name);
    }

    auto(info: FileInfo, dataStateObject: PluginStateObject.Data.Binary | PluginStateObject.Data.String) {
        for (let i = 0, il = this.list.length; i < il; ++i) {
            const { provider } = this._list[i];

            let hasExt = false;
            if (provider.binaryExtensions && provider.binaryExtensions.indexOf(info.ext) >= 0) hasExt = true;
            else if (provider.stringExtensions && provider.stringExtensions.indexOf(info.ext) >= 0) hasExt = true;
            if (hasExt && (!provider.isApplicable || provider.isApplicable(info, dataStateObject.data))) return provider;
        }
        return;
    }

    get(name: string): DataFormatProvider | undefined {
        if (this._map.has(name)) {
            return this._map.get(name)!;
        } else {
            throw new Error(`unknown data format name '${name}'`);
        }
    }

    get list() {
        return this._list;
    }
}