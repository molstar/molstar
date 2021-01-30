/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import '../../mol-util/polyfill';
import { createPlugin, DefaultPluginSpec } from '../../mol-plugin';
import './index.html';
import './embedded.html';
import './favicon.ico';
import { PluginContext } from '../../mol-plugin/context';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginSpec } from '../../mol-plugin/spec';
import { DownloadStructure, PdbDownloadProvider } from '../../mol-plugin-state/actions/structure';
import { PluginConfig } from '../../mol-plugin/config';
import { CellPack } from '../../extensions/cellpack';
import { RCSBAssemblySymmetry, RCSBValidationReport } from '../../extensions/rcsb';
import { PDBeStructureQualityReport } from '../../extensions/pdbe';
import { Asset } from '../../mol-util/assets';
import { ObjectKeys } from '../../mol-util/type-helpers';
import { PluginState } from '../../mol-plugin/state';
import { DownloadDensity } from '../../mol-plugin-state/actions/volume';
import { PluginLayoutControlsDisplay } from '../../mol-plugin/layout';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { ANVILMembraneOrientation } from '../../extensions/anvil/behavior';
import { DnatcoConfalPyramids } from '../../extensions/dnatco';
import { G3DFormat, G3dProvider } from '../../extensions/g3d/format';
import { DataFormatProvider } from '../../mol-plugin-state/formats/provider';
import { BuildInVolumeFormat } from '../../mol-plugin-state/formats/volume';
import { Color } from '../../mol-util/color';
import { StateObjectSelector } from '../../mol-state';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { Mp4Export } from '../../extensions/mp4-export';
import { StructureRepresentationPresetProvider } from '../../mol-plugin-state/builder/structure/representation-preset';

require('mol-plugin-ui/skin/light.scss');

export { PLUGIN_VERSION as version } from '../../mol-plugin/version';
export { setProductionMode, setDebugMode } from '../../mol-util/debug';

const CustomFormats = [
    ['g3d', G3dProvider] as const
];

const Extensions = {
    'cellpack': PluginSpec.Behavior(CellPack),
    'dnatco-confal-pyramids': PluginSpec.Behavior(DnatcoConfalPyramids),
    'pdbe-structure-quality-report': PluginSpec.Behavior(PDBeStructureQualityReport),
    'rcsb-assembly-symmetry': PluginSpec.Behavior(RCSBAssemblySymmetry),
    'rcsb-validation-report': PluginSpec.Behavior(RCSBValidationReport),
    'anvil-membrane-orientation': PluginSpec.Behavior(ANVILMembraneOrientation),
    'g3d': PluginSpec.Behavior(G3DFormat),
    'mp4-export': PluginSpec.Behavior(Mp4Export)
};

const DefaultViewerOptions = {
    customFormats: CustomFormats as [string, DataFormatProvider][],
    extensions: ObjectKeys(Extensions),
    layoutIsExpanded: true,
    layoutShowControls: true,
    layoutShowRemoteState: true,
    layoutControlsDisplay: 'reactive' as PluginLayoutControlsDisplay,
    layoutShowSequence: true,
    layoutShowLog: true,
    layoutShowLeftPanel: true,
    disableAntialiasing: false,
    pixelScale: 1,
    enableWboit: true,

    viewportShowExpand: PluginConfig.Viewport.ShowExpand.defaultValue,
    viewportShowControls: PluginConfig.Viewport.ShowControls.defaultValue,
    viewportShowSettings: PluginConfig.Viewport.ShowSettings.defaultValue,
    viewportShowSelectionMode: PluginConfig.Viewport.ShowSelectionMode.defaultValue,
    viewportShowAnimation: PluginConfig.Viewport.ShowAnimation.defaultValue,
    pluginStateServer: PluginConfig.State.DefaultServer.defaultValue,
    volumeStreamingServer: PluginConfig.VolumeStreaming.DefaultServer.defaultValue,
    volumeStreamingDisabled: !PluginConfig.VolumeStreaming.Enabled.defaultValue,
    pdbProvider: PluginConfig.Download.DefaultPdbProvider.defaultValue,
    emdbProvider: PluginConfig.Download.DefaultEmdbProvider.defaultValue,
};
type ViewerOptions = typeof DefaultViewerOptions;

export class Viewer {
    plugin: PluginContext

    constructor(elementOrId: string | HTMLElement, options: Partial<ViewerOptions> = {}) {
        const o = { ...DefaultViewerOptions, ...options };

        const spec: PluginSpec = {
            actions: [...DefaultPluginSpec.actions],
            behaviors: [
                ...DefaultPluginSpec.behaviors,
                ...o.extensions.map(e => Extensions[e]),
            ],
            animations: [...DefaultPluginSpec.animations || []],
            customParamEditors: DefaultPluginSpec.customParamEditors,
            customFormats: o?.customFormats,
            layout: {
                initial: {
                    isExpanded: o.layoutIsExpanded,
                    showControls: o.layoutShowControls,
                    controlsDisplay: o.layoutControlsDisplay,
                },
                controls: {
                    ...DefaultPluginSpec.layout && DefaultPluginSpec.layout.controls,
                    top: o.layoutShowSequence ? undefined : 'none',
                    bottom: o.layoutShowLog ? undefined : 'none',
                    left: o.layoutShowLeftPanel ? undefined : 'none',
                }
            },
            components: {
                ...DefaultPluginSpec.components,
                remoteState: o.layoutShowRemoteState ? 'default' : 'none',
            },
            config: [
                [PluginConfig.General.DisableAntialiasing, o.disableAntialiasing],
                [PluginConfig.General.PixelScale, o.pixelScale],
                [PluginConfig.General.EnableWboit, o.enableWboit],
                [PluginConfig.Viewport.ShowExpand, o.viewportShowExpand],
                [PluginConfig.Viewport.ShowControls, o.viewportShowControls],
                [PluginConfig.Viewport.ShowSettings, o.viewportShowSettings],
                [PluginConfig.Viewport.ShowSelectionMode, o.viewportShowSelectionMode],
                [PluginConfig.Viewport.ShowAnimation, o.viewportShowAnimation],
                [PluginConfig.State.DefaultServer, o.pluginStateServer],
                [PluginConfig.State.CurrentServer, o.pluginStateServer],
                [PluginConfig.VolumeStreaming.DefaultServer, o.volumeStreamingServer],
                [PluginConfig.VolumeStreaming.Enabled, !o.volumeStreamingDisabled],
                [PluginConfig.Download.DefaultPdbProvider, o.pdbProvider],
                [PluginConfig.Download.DefaultEmdbProvider, o.emdbProvider]
            ]
        };

        const element = typeof elementOrId === 'string'
            ? document.getElementById(elementOrId)
            : elementOrId;
        if (!element) throw new Error(`Could not get element with id '${elementOrId}'`);
        this.plugin = createPlugin(element, spec);
    }

    setRemoteSnapshot(id: string) {
        const url = `${this.plugin.config.get(PluginConfig.State.CurrentServer)}/get/${id}`;
        return PluginCommands.State.Snapshots.Fetch(this.plugin, { url });
    }

    loadSnapshotFromUrl(url: string, type: PluginState.SnapshotType) {
        return PluginCommands.State.Snapshots.OpenUrl(this.plugin, { url, type });
    }

    loadStructureFromUrl(url: string, format: BuiltInTrajectoryFormat = 'mmcif', isBinary = false, options?: LoadStructureOptions) {
        const params = DownloadStructure.createDefaultParams(this.plugin.state.data.root.obj!, this.plugin);
        return this.plugin.runTask(this.plugin.state.data.applyAction(DownloadStructure, {
            source: {
                name: 'url',
                params: {
                    url: Asset.Url(url),
                    format: format as any,
                    isBinary,
                    options: { ...params.source.params.options, representationParams: options?.representationParams as any },
                }
            }
        }));
    }

    async loadAllModelsOrAssemblyFromUrl(url: string, format: BuiltInTrajectoryFormat = 'mmcif', isBinary = false, options?: LoadStructureOptions) {
        const plugin = this.plugin;

        const data = await plugin.builders.data.download({ url, isBinary }, { state: { isGhost: true } });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, format);

        await this.plugin.builders.structure.hierarchy.applyPreset(trajectory, 'all-models', { useDefaultIfSingleModel: true, representationPresetParams: options?.representationParams });
    }

    async loadStructureFromData(data: string | number[], format: BuiltInTrajectoryFormat, options?: { dataLabel?: string }) {
        const _data = await this.plugin.builders.data.rawData({ data, label: options?.dataLabel });
        const trajectory = await this.plugin.builders.structure.parseTrajectory(_data, format);
        await this.plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    }

    loadPdb(pdb: string, options?: LoadStructureOptions) {
        const params = DownloadStructure.createDefaultParams(this.plugin.state.data.root.obj!, this.plugin);
        const provider = this.plugin.config.get(PluginConfig.Download.DefaultPdbProvider)!;
        return this.plugin.runTask(this.plugin.state.data.applyAction(DownloadStructure, {
            source: {
                name: 'pdb' as const,
                params: {
                    provider: {
                        id: pdb,
                        server: {
                            name: provider,
                            params: PdbDownloadProvider[provider].defaultValue as any
                        }
                    },
                    options: { ...params.source.params.options, representationParams: options?.representationParams as any },
                }
            }
        }));
    }

    loadPdbDev(pdbDev: string) {
        const params = DownloadStructure.createDefaultParams(this.plugin.state.data.root.obj!, this.plugin);
        return this.plugin.runTask(this.plugin.state.data.applyAction(DownloadStructure, {
            source: {
                name: 'pdb-dev' as const,
                params: {
                    provider: {
                        id: pdbDev,
                        encoding: 'bcif',
                    },
                    options: params.source.params.options,
                }
            }
        }));
    }

    loadEmdb(emdb: string, options?: { detail?: number }) {
        const provider = this.plugin.config.get(PluginConfig.Download.DefaultEmdbProvider)!;
        return this.plugin.runTask(this.plugin.state.data.applyAction(DownloadDensity, {
            source: {
                name: 'pdb-emd-ds' as const,
                params: {
                    provider: {
                        id: emdb,
                        server: provider,
                    },
                    detail: options?.detail ?? 3,
                }
            }
        }));
    }

    async loadVolumeFromUrl(url: string, format: BuildInVolumeFormat, isBinary: boolean, isovalues: VolumeIsovalueInfo[], entryId?: string) {
        const plugin = this.plugin;

        if (!plugin.dataFormats.get(format)) {
            throw new Error(`Unknown density format: ${format}`);
        }

        return plugin.dataTransaction(async () => {
            const data = await plugin.builders.data.download({ url, isBinary, label: entryId }, { state: { isGhost: true } });

            const parsed = await plugin.dataFormats.get(format)!.parse(plugin, data, { entryId });
            const volume = (parsed.volume || parsed.volumes[0]) as StateObjectSelector<PluginStateObject.Volume.Data>;
            if (!volume?.isOk) throw new Error('Failed to parse any volume.');

            const repr = plugin.build().to(volume);
            for (const iso of isovalues) {
                repr.apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.plugin, volume.data!, {
                    type: 'isosurface',
                    typeParams: { alpha: iso.alpha ?? 1, isoValue: iso.type === 'absolute' ?  { kind: 'absolute', absoluteValue: iso.value } : { kind: 'relative', relativeValue: iso.value } },
                    color: 'uniform',
                    colorParams: { value: iso.color }
                }));
            }

            await repr.commit();
        });
    }

    handleResize() {
        this.plugin.layout.events.updated.next();
    }
}

export interface LoadStructureOptions {
    representationParams?: StructureRepresentationPresetProvider.CommonParams
}

export interface VolumeIsovalueInfo {
    type: 'absolute' | 'relative',
    value: number,
    color: Color,
    alpha?: number
}