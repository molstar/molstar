/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export class PluginConfigItem<T = any> {
    toString() { return this.key; }
    valueOf() { return this.key; }
    constructor(public key: string, public defaultValue?: T) { }
}

function item<T>(key: string, defaultValue?: T) { return new PluginConfigItem(key, defaultValue); }

export const PluginConfig = {
    item,
    PluginState: { Server: item('plugin-state.server', 'https://webchem.ncbr.muni.cz/molstar-state') }
}

export class PluginConfigManager {    
    private _config = new Map<PluginConfigItem<any>, unknown>();

    get<T>(key: PluginConfigItem<T>) {
        if (!this._config.has(key)) return key.defaultValue;
        return this._config.get(key) as T;
    }

    set<T>(key: PluginConfigItem<T>, value: T) {
        this._config.set(key, value);
    }

    delete<T>(key: PluginConfigItem<T>) {
        this._config.delete(key);
    }

    constructor(initial?: Map<PluginConfigItem, unknown>) {
        if (!initial) return;
        initial.forEach((v, k) => this._config.set(k, v));
    }
}