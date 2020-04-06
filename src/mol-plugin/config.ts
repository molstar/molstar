/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Model } from '../mol-model/structure';
import { PluginContext } from './context';

export class PluginConfigItem<T = any> {
    toString() { return this.key; }
    valueOf() { return this.key; }
    constructor(public key: string, public defaultValue?: T) { }
}

function item<T>(key: string, defaultValue?: T) { return new PluginConfigItem(key, defaultValue); }

export const PluginConfig = {
    item,
    General: {
        IsBusyTimeoutMs: item('plugin-config.is-busy-timeout', 750)
    },
    State: {
        DefaultServer: item('plugin-state.server', 'https://webchem.ncbr.muni.cz/molstar-state'),
        CurrentServer: item('plugin-state.server', 'https://webchem.ncbr.muni.cz/molstar-state')
    },
    VolumeStreaming: {
        DefaultServer: item('volume-streaming.server', 'https://ds.litemol.org'),
        CanStream: item('volume-streaming.can-stream', (s: Structure, plugin: PluginContext) => {
            return s.models.length === 1 && (Model.hasDensityMap(s.models[0])
                // the following test is to include e.g. 'updated' files from PDBe
                || (!Model.isFromPdbArchive(s.models[0]) && s.models[0].entryId.length === 4))
        }),
        EmdbHeaderServer: item('volume-streaming.emdb-header-server', 'https://ftp.wwpdb.org/pub/emdb/structures'),
    },
    Viewport: {
        ShowExpand: item('viewer.show-expand-button', true),
        ShowSelectionMode: item('viewer.show-selection-model-button', true)
    }
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