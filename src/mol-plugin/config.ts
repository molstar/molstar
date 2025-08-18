/**
 * Copyright (c) 2020-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Model } from '../mol-model/structure';
import { PluginContext } from './context';
import { PdbDownloadProvider } from '../mol-plugin-state/actions/structure';
import { EmdbDownloadProvider } from '../mol-plugin-state/actions/volume';
import { StructureRepresentationPresetProvider } from '../mol-plugin-state/builder/structure/representation-preset';
import { PluginFeatureDetection } from './features';
import { SaccharideCompIdMapType } from '../mol-model/structure/structure/carbohydrates/constants';
import { BackgroundProps } from '../mol-canvas3d/passes/background';

export class PluginConfigItem<T = any> {
    toString() { return this.key; }
    valueOf() { return this.key; }
    constructor(public key: string, public defaultValue?: T) { }
}

function item<T>(key: string, defaultValue?: T) { return new PluginConfigItem(key, defaultValue); }

export const PluginConfig = {
    item,
    General: {
        IsBusyTimeoutMs: item('plugin-config.is-busy-timeout', 750),
        DisableAntialiasing: item('plugin-config.disable-antialiasing', false),
        DisablePreserveDrawingBuffer: item('plugin-config.disable-preserve-drawing-buffer', false),
        PixelScale: item('plugin-config.pixel-scale', 1),
        PickScale: item('plugin-config.pick-scale', 0.25),
        Transparency: item<'blended' | 'wboit' | 'dpoit'>('plugin-config.transparency', PluginFeatureDetection.defaultTransparency),
        // as of Oct 1 2021, WebGL 2 doesn't work on iOS 15.
        // TODO: check back in a few weeks to see if it was fixed
        PreferWebGl1: item('plugin-config.prefer-webgl1', PluginFeatureDetection.preferWebGl1),
        AllowMajorPerformanceCaveat: item('plugin-config.allow-major-performance-caveat', false),
        PowerPreference: item<WebGLContextAttributes['powerPreference']>('plugin-config.power-preference', 'high-performance'),
        ResolutionMode: item<'auto' | 'scaled' | 'native'>('plugin-config.resolution-mode', 'auto'),
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
        EmdbHeaderServer: item('volume-streaming.emdb-header-server', 'https://files.wwpdb.org/pub/emdb/structures'),
    },
    Viewport: {
        ShowReset: item('viewer.show-reset-button', true),
        ShowExpand: item('viewer.show-expand-button', true),
        ShowControls: item('viewer.show-controls-button', true),
        ShowSettings: item('viewer.show-settings-button', true),
        ShowSelectionMode: item('viewer.show-selection-model-button', true),
        ShowAnimation: item('viewer.show-animation-button', true),
        ShowTrajectoryControls: item('viewer.show-trajectory-controls', true),
        ShowScreenshotControls: item('viewer.show-screenshot-controls', true),
    },
    Download: {
        DefaultPdbProvider: item<PdbDownloadProvider>('download.default-pdb-provider', 'pdbe'),
        DefaultEmdbProvider: item<EmdbDownloadProvider>('download.default-emdb-provider', 'pdbe'),
    },
    Structure: {
        SizeThresholds: item('structure.size-thresholds', Structure.DefaultSizeThresholds),
        DefaultRepresentationPreset: item<string>('structure.default-representation-preset', 'auto'),
        DefaultRepresentationPresetParams: item<StructureRepresentationPresetProvider.CommonParams>('structure.default-representation-preset-params', { }),
        SaccharideCompIdMapType: item<SaccharideCompIdMapType>('structure.saccharide-comp-id-map-type', 'default'),
    },
    Background: {
        Styles: item<[BackgroundProps, string][]>('background.styles', []),
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
