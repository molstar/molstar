/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Model } from '../mol-model/structure';
import { PluginContext } from './context';
import { PdbDownloadProvider } from '../mol-plugin-state/actions/structure';
import { EmdbDownloadProvider } from '../mol-plugin-state/actions/volume';
import { StructureRepresentationPresetProvider } from '../mol-plugin-state/builder/structure/representation-preset';

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
            return s.models.length === 1 && Model.probablyHasDensityMap(s.models[0]);
        }),
        EmdbHeaderServer: item('volume-streaming.emdb-header-server', 'https://ftp.wwpdb.org/pub/emdb/structures'),
    },
    Viewport: {
        ShowExpand: item('viewer.show-expand-button', true),
        ShowSelectionMode: item('viewer.show-selection-model-button', true),
        ShowAnimation: item('viewer.show-animation-button', true),
    },
    Download: {
        DefaultPdbProvider: item<PdbDownloadProvider>('download.default-pdb-provider', 'pdbe'),
        DefaultEmdbProvider: item<EmdbDownloadProvider>('download.default-emdb-provider', 'pdbe'),
    },
    Structure: {
        SizeThresholds: item('structure.size-thresholds', Structure.DefaultSizeThresholds),
        DefaultRepresentationPresetParams: item<StructureRepresentationPresetProvider.CommonParams>('structure.default-representation-preset-params', { })
    }
};

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

    constructor(initial?: [PluginConfigItem, unknown][]) {
        if (!initial) return;
        initial.forEach(([k, v]) => this._config.set(k, v));
    }
}