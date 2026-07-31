/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FileNameInfo } from '../../mol-util/file-info';
import { PluginStateObject } from '../objects';
import { DataFormatProvider } from './provider';
import { BuiltInTrajectoryFormats } from './trajectory';
import { BuiltInVolumeFormats } from './volume';
import { BuiltInShapeFormats } from './shape';
import { BuiltInTopologyFormats } from './topology';
import { BuiltInCoordinatesFormats } from './coordinates';

export class DataFormatRegistry {
    private _list: { name: string, provider: DataFormatProvider }[] = [];
    private _map = new Map<string, DataFormatProvider>();
    private _extensions: Set<string> | undefined = undefined;
    private _binaryExtensions: Set<string> | undefined = undefined;
    private _options: [name: string, label: string, category: string][] | undefined = undefined;
    private _autoOrder: { name: string, provider: DataFormatProvider }[] | undefined = undefined;

    get types(): [name: string, label: string][] {
        return this._list.map(e => [e.name, e.provider.label] as [name: string, label: string]);
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
        const options: [name: string, label: string, category: string][] = [];
        this._list.forEach(({ name, provider }) => options.push([name, provider.label, provider.category || '']));
        this._options = options;
        return options;
    }

    /**
     * Providers ordered for `auto()`: higher `priority` first, ties broken by registration
     * order (stable sort). Explicit priority makes auto-detection independent of the order
     * providers happen to be registered in the constructor below.
     */
    get autoOrder() {
        if (this._autoOrder) return this._autoOrder;
        this._autoOrder = this._list
            .map((entry, index) => ({ entry, index }))
            .sort((a, b) => (b.entry.provider.priority ?? 0) - (a.entry.provider.priority ?? 0) || a.index - b.index)
            .map(({ entry }) => entry);
        return this._autoOrder;
    }

    constructor() {
        for (const [id, p] of BuiltInVolumeFormats) this.add(id, p);
        for (const [id, p] of BuiltInTopologyFormats) this.add(id, p);
        for (const [id, p] of BuiltInCoordinatesFormats) this.add(id, p);
        for (const [id, p] of BuiltInShapeFormats) this.add(id, p);
        for (const [id, p] of BuiltInTrajectoryFormats) this.add(id, p);
    };

    private _clear() {
        this._extensions = undefined;
        this._binaryExtensions = undefined;
        this._options = undefined;
        this._autoOrder = undefined;
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

    auto(info: FileNameInfo, dataStateObject: PluginStateObject.Data.Binary | PluginStateObject.Data.String) {
        const list = this.autoOrder;
        for (let i = 0, il = list.length; i < il; ++i) {
            const p = list[i].provider;

            let hasExt = false;
            if (p.binaryExtensions?.includes(info.ext)) hasExt = true;
            else if (p.stringExtensions?.includes(info.ext)) hasExt = true;
            if (hasExt && (!p.isApplicable || p.isApplicable(info, dataStateObject.data))) return p;
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