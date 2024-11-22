/**
 * Copyright (c) 2022-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mp4Export } from '../../extensions/mp4-export';
import { DataFormatProvider } from '../../mol-plugin-state/formats/provider';
import { createPluginUI } from '../../mol-plugin-ui';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { PluginUIContext } from '../../mol-plugin-ui/context';
import { DefaultPluginUISpec, PluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginLayoutControlsDisplay } from '../../mol-plugin/layout';
import { PluginSpec } from '../../mol-plugin/spec';
import '../../mol-util/polyfill';
import { ObjectKeys } from '../../mol-util/type-helpers';
import { SaccharideCompIdMapType } from '../../mol-model/structure/structure/carbohydrates/constants';
import { Backgrounds } from '../../extensions/backgrounds';
import { LeftPanel, RightPanel } from './ui/panels';
import { Color } from '../../mol-util/color';
import { SpacefillRepresentationProvider } from '../../mol-repr/structure/representation/spacefill';
import { PluginBehaviors } from '../../mol-plugin/behavior';
import { MesoFocusLoci } from './behavior/camera';
import { GraphicsMode, MesoscaleState } from './data/state';
import { MesoSelectLoci } from './behavior/select';
import { Transparency } from '../../mol-gl/webgl/render-item';
import { LoadModel, loadExampleEntry, loadPdb, loadPdbDev, loadUrl, openState } from './ui/states';
import { Asset } from '../../mol-util/assets';
import { AnimateCameraSpin } from '../../mol-plugin-state/animation/built-in/camera-spin';
import { AnimateCameraRock } from '../../mol-plugin-state/animation/built-in/camera-rock';
import { AnimateStateSnapshots } from '../../mol-plugin-state/animation/built-in/state-snapshots';
import { MesoViewportSnapshotDescription } from './ui/entities';

export { PLUGIN_VERSION as version } from '../../mol-plugin/version';
export { setDebugMode, setProductionMode, setTimingMode, consoleStats } from '../../mol-util/debug';

export type ExampleEntry = {
    id: string,
    label: string,
    url: string,
    type: 'molx' | 'molj' | 'cif' | 'bcif',
    description?: string,
    link?: string,
}

export type MesoscaleExplorerState = {
    examples?: ExampleEntry[],
    graphicsMode: GraphicsMode,
    illumination: boolean,
    stateRef?: string,
    driver?: any,
    stateCache: { [k: string]: any },
}

//

const Extensions = {
    'backgrounds': PluginSpec.Behavior(Backgrounds),
    'mp4-export': PluginSpec.Behavior(Mp4Export),
};

const DefaultMesoscaleExplorerOptions = {
    customFormats: [] as [string, DataFormatProvider][],
    extensions: ObjectKeys(Extensions),
    layoutIsExpanded: true,
    layoutShowControls: true,
    layoutShowRemoteState: true,
    layoutControlsDisplay: 'reactive' as PluginLayoutControlsDisplay,
    layoutShowSequence: true,
    layoutShowLog: true,
    layoutShowLeftPanel: true,
    collapseLeftPanel: false,
    collapseRightPanel: false,
    disableAntialiasing: PluginConfig.General.DisableAntialiasing.defaultValue,
    pixelScale: PluginConfig.General.PixelScale.defaultValue,
    pickScale: PluginConfig.General.PickScale.defaultValue,
    transparency: 'blended' as Transparency,
    preferWebgl1: PluginConfig.General.PreferWebGl1.defaultValue,
    allowMajorPerformanceCaveat: PluginConfig.General.AllowMajorPerformanceCaveat.defaultValue,
    powerPreference: PluginConfig.General.PowerPreference.defaultValue,
    resolutionMode: PluginConfig.General.ResolutionMode.defaultValue,
    illumination: false,

    viewportShowExpand: PluginConfig.Viewport.ShowExpand.defaultValue,
    viewportShowControls: PluginConfig.Viewport.ShowControls.defaultValue,
    viewportShowSettings: PluginConfig.Viewport.ShowSettings.defaultValue,
    viewportShowSelectionMode: false,
    viewportShowAnimation: false,
    viewportShowTrajectoryControls: false,
    pluginStateServer: PluginConfig.State.DefaultServer.defaultValue,
    volumeStreamingServer: PluginConfig.VolumeStreaming.DefaultServer.defaultValue,
    volumeStreamingDisabled: !PluginConfig.VolumeStreaming.Enabled.defaultValue,
    pdbProvider: PluginConfig.Download.DefaultPdbProvider.defaultValue,
    emdbProvider: PluginConfig.Download.DefaultEmdbProvider.defaultValue,
    saccharideCompIdMapType: 'default' as SaccharideCompIdMapType,

    graphicsMode: 'quality' as GraphicsMode,
    driver: undefined
};
type MesoscaleExplorerOptions = typeof DefaultMesoscaleExplorerOptions;

export class MesoscaleExplorer {
    constructor(public plugin: PluginUIContext) {
    }

    async loadExample(id: string) {
        const entries = (this.plugin.customState as MesoscaleExplorerState).examples || [];
        const entry = entries.find(e => e.id === id);
        if (entry !== undefined) {
            await loadExampleEntry(this.plugin, entry);
        }
    }

    async loadUrl(url: string, type: 'molx' | 'molj' | 'cif' | 'bcif') {
        await loadUrl(this.plugin, url, type);
    }

    async loadPdb(id: string) {
        await loadPdb(this.plugin, id);
    }

    async loadPdbDev(id: string) {
        await loadPdbDev(this.plugin, id);
    }

    static async create(elementOrId: string | HTMLElement, options: Partial<MesoscaleExplorerOptions> = {}) {
        const definedOptions = {} as any;
        // filter for defined properies only so the default values
        // are property applied
        for (const p of Object.keys(options) as (keyof MesoscaleExplorerOptions)[]) {
            if (options[p] !== void 0) definedOptions[p] = options[p];
        }

        const o: MesoscaleExplorerOptions = { ...DefaultMesoscaleExplorerOptions, ...definedOptions };
        const defaultSpec = DefaultPluginUISpec();

        const spec: PluginUISpec = {
            actions: defaultSpec.actions,
            behaviors: [
                PluginSpec.Behavior(PluginBehaviors.Camera.CameraAxisHelper),
                PluginSpec.Behavior(PluginBehaviors.Camera.CameraControls),

                PluginSpec.Behavior(MesoFocusLoci),
                PluginSpec.Behavior(MesoSelectLoci),

                ...o.extensions.map(e => Extensions[e]),
            ],
            animations: [
                AnimateCameraSpin,
                AnimateCameraRock,
                AnimateStateSnapshots,
            ],
            customParamEditors: defaultSpec.customParamEditors,
            customFormats: o?.customFormats,
            layout: {
                initial: {
                    isExpanded: o.layoutIsExpanded,
                    showControls: o.layoutShowControls,
                    controlsDisplay: o.layoutControlsDisplay,
                    regionState: {
                        bottom: 'full',
                        left: o.collapseLeftPanel ? 'collapsed' : 'full',
                        right: o.collapseRightPanel ? 'hidden' : 'full',
                        top: 'full',
                    }
                },
            },
            components: {
                ...defaultSpec.components,
                controls: {
                    ...defaultSpec.components?.controls,
                    top: 'none',
                    bottom: 'none',
                    left: LeftPanel,
                    right: RightPanel,
                },
                remoteState: 'none',
                viewport: {
                    snapshotDescription: MesoViewportSnapshotDescription,
                }
            },
            config: [
                [PluginConfig.General.DisableAntialiasing, o.disableAntialiasing],
                [PluginConfig.General.PixelScale, o.pixelScale],
                [PluginConfig.General.PickScale, o.pickScale],
                [PluginConfig.General.Transparency, o.transparency],
                [PluginConfig.General.PreferWebGl1, o.preferWebgl1],
                [PluginConfig.General.AllowMajorPerformanceCaveat, o.allowMajorPerformanceCaveat],
                [PluginConfig.General.PowerPreference, o.powerPreference],
                [PluginConfig.General.ResolutionMode, o.resolutionMode],
                [PluginConfig.Viewport.ShowExpand, o.viewportShowExpand],
                [PluginConfig.Viewport.ShowControls, o.viewportShowControls],
                [PluginConfig.Viewport.ShowSettings, o.viewportShowSettings],
                [PluginConfig.Viewport.ShowSelectionMode, o.viewportShowSelectionMode],
                [PluginConfig.Viewport.ShowAnimation, o.viewportShowAnimation],
                [PluginConfig.Viewport.ShowTrajectoryControls, o.viewportShowTrajectoryControls],
                [PluginConfig.State.DefaultServer, o.pluginStateServer],
                [PluginConfig.State.CurrentServer, o.pluginStateServer],
                [PluginConfig.VolumeStreaming.DefaultServer, o.volumeStreamingServer],
                [PluginConfig.VolumeStreaming.Enabled, !o.volumeStreamingDisabled],
                [PluginConfig.Download.DefaultPdbProvider, o.pdbProvider],
                [PluginConfig.Download.DefaultEmdbProvider, o.emdbProvider],
                [PluginConfig.Structure.SaccharideCompIdMapType, o.saccharideCompIdMapType],
            ]
        };

        const element = typeof elementOrId === 'string'
            ? document.getElementById(elementOrId)
            : elementOrId;
        if (!element) throw new Error(`Could not get element with id '${elementOrId}'`);

        const plugin = await createPluginUI({
            target: element,
            spec,
            render: renderReact18,
            onBeforeUIRender: async plugin => {
                let examples: MesoscaleExplorerState['examples'] = undefined;
                try {
                    examples = await plugin.fetch({ url: '../examples/list.json', type: 'json' }).run();
                    // extend the array with file tour.json if it exists
                    const tour = await plugin.fetch({ url: '../examples/tour.json', type: 'json' }).run();
                    if (tour) {
                        examples = examples?.concat(tour);
                    }
                } catch (e) {
                    console.log(e);
                }

                (plugin.customState as MesoscaleExplorerState) = {
                    examples,
                    graphicsMode: o.graphicsMode,
                    illumination: o.illumination,
                    driver: o.driver,
                    stateCache: {},
                };

                await MesoscaleState.init(plugin);
            }
        });

        plugin.canvas3d?.setProps({
            renderer: {
                backgroundColor: Color(0x101010),
            },
            cameraFog: { name: 'off', params: {} },
            hiZ: { enabled: true },
        });

        plugin.representation.structure.registry.clear();
        plugin.representation.structure.registry.add(SpacefillRepresentationProvider);

        plugin.state.setSnapshotParams({
            image: true,
            componentManager: false,
            structureSelection: true,
            behavior: true,
        });

        plugin.managers.lociLabels.clearProviders();

        plugin.managers.dragAndDrop.addHandler('mesoscale-explorer', (files) => {
            const sessions = files.filter(f => {
                const fn = f.name.toLowerCase();
                return fn.endsWith('.molx') || fn.endsWith('.molj');
            });

            if (sessions.length > 0) {
                openState(plugin, sessions[0]);
            } else {
                plugin.runTask(plugin.state.data.applyAction(LoadModel, {
                    files: files.map(f => Asset.File(f)),
                }));
            }

            return true;
        });

        plugin.state.events.object.created.subscribe(e => {
            (plugin.customState as MesoscaleExplorerState).stateCache = {};
        });

        plugin.state.events.object.removed.subscribe(e => {
            (plugin.customState as MesoscaleExplorerState).stateCache = {};
        });

        return new MesoscaleExplorer(plugin);
    }

    handleResize() {
        this.plugin.layout.events.updated.next(void 0);
    }

    dispose() {
        this.plugin.dispose();
    }
}
