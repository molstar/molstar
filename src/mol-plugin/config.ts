/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
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


function preferWebGl1() {
    if (typeof navigator === 'undefined' || typeof window === 'undefined') return false;

    // WebGL2 isn't working in MacOS 12.0.1 Safari 15.1 (but is working in Safari tech preview)
    // prefer webgl 1 based on the userAgent substring
    if (navigator.userAgent.indexOf('Version/15.1 Safari') > 0) {
        return true;
    }

    // Check for iOS device which enabled WebGL2 recently but it doesn't seem
    // to be full up to speed yet.

    // adapted from https://stackoverflow.com/questions/9038625/detect-if-device-is-ios
    const isIOS = /iPad|iPhone|iPod/.test(navigator.userAgent);
    const isAppleDevice = navigator.userAgent.includes('Macintosh');
    const isTouchScreen = navigator.maxTouchPoints >= 4; // true for iOS 13 (and hopefully beyond)
    return !(window as any).MSStream && (isIOS || (isAppleDevice && isTouchScreen));
}

export const PluginConfig = {
    item,
    General: {
        IsBusyTimeoutMs: item('plugin-config.is-busy-timeout', 750),
        DisableAntialiasing: item('plugin-config.disable-antialiasing', false),
        DisablePreserveDrawingBuffer: item('plugin-config.disable-preserve-drawing-buffer', false),
        PixelScale: item('plugin-config.pixel-scale', 1),
        PickScale: item('plugin-config.pick-scale', 0.25),
        PickPadding: item('plugin-config.pick-padding', 3),
        EnableWboit: item('plugin-config.enable-wboit', true),
        // as of Oct 1 2021, WebGL 2 doesn't work on iOS 15.
        // TODO: check back in a few weeks to see if it was fixed
        PreferWebGl1: item('plugin-config.prefer-webgl1', preferWebGl1()),
    },
    State: {
        DefaultServer: item('plugin-state.server', 'https://webchem.ncbr.muni.cz/molstar-state'),
        CurrentServer: item('plugin-state.server', 'https://webchem.ncbr.muni.cz/molstar-state'),
        HistoryCapacity: item('history-capacity.server', 5)
    },
    VolumeStreaming: {
        Enabled: item('volume-streaming.enabled', true),
        DefaultServer: item('volume-streaming.server', 'https://ds.litemol.org'),
        CanStream: item('volume-streaming.can-stream', (s: Structure, plugin: PluginContext) => {
            return s.models.length === 1 && Model.probablyHasDensityMap(s.models[0]);
        }),
        EmdbHeaderServer: item('volume-streaming.emdb-header-server', 'https://ftp.wwpdb.org/pub/emdb/structures'),
    },
    Viewport: {
        ShowExpand: item('viewer.show-expand-button', true),
        ShowControls: item('viewer.show-controls-button', true),
        ShowSettings: item('viewer.show-settings-button', true),
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