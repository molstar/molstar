/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import '../../mol-util/polyfill';
import { createPlugin, DefaultPluginSpec } from '../../mol-plugin';
import './index.html';
import { PluginContext } from '../../mol-plugin/context';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginSpec } from '../../mol-plugin/spec';
import { DownloadStructure, PdbDownloadProvider } from '../../mol-plugin-state/actions/structure';
import { PluginConfig } from '../../mol-plugin/config';
import { Asset } from '../../mol-util/assets';
import { ObjectKeys } from '../../mol-util/type-helpers';
import { PluginState } from '../../mol-plugin/state';
import { DownloadDensity } from '../../mol-plugin-state/actions/volume';
import { PluginLayoutControlsDisplay } from '../../mol-plugin/layout';
import { BuiltInTrajectoryFormat } from '../../mol-plugin-state/formats/trajectory';
import { Structure } from '../../mol-model/structure';
import { PluginStateTransform, PluginStateObject as PSO } from '../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Task } from '../../mol-task';
import { StateObject } from '../../mol-state';
import { ViewportComponent, StructurePreset } from './viewport';
import { PluginBehaviors } from '../../mol-plugin/behavior';
import { ColorNames } from '../../mol-util/color/names';

require('mol-plugin-ui/skin/light.scss');

export { PLUGIN_VERSION as version } from '../../mol-plugin/version';
export { setProductionMode, setDebugMode } from '../../mol-util/debug';

const DefaultViewerOptions = {
    extensions: ObjectKeys({}),
    layoutIsExpanded: true,
    layoutShowControls: true,
    layoutShowRemoteState: true,
    layoutControlsDisplay: 'reactive' as PluginLayoutControlsDisplay,
    layoutShowSequence: true,
    layoutShowLog: true,
    layoutShowLeftPanel: true,

    viewportShowExpand: PluginConfig.Viewport.ShowExpand.defaultValue,
    viewportShowControls: PluginConfig.Viewport.ShowControls.defaultValue,
    viewportShowSettings: PluginConfig.Viewport.ShowSettings.defaultValue,
    viewportShowSelectionMode: PluginConfig.Viewport.ShowSelectionMode.defaultValue,
    viewportShowAnimation: PluginConfig.Viewport.ShowAnimation.defaultValue,
    pluginStateServer: PluginConfig.State.DefaultServer.defaultValue,
    volumeStreamingServer: PluginConfig.VolumeStreaming.DefaultServer.defaultValue,
    pdbProvider: PluginConfig.Download.DefaultPdbProvider.defaultValue,
    emdbProvider: PluginConfig.Download.DefaultEmdbProvider.defaultValue,
};
type ViewerOptions = typeof DefaultViewerOptions;

class Viewer {
    plugin: PluginContext

    constructor(elementOrId: string | HTMLElement, options: Partial<ViewerOptions> = {}) {
        const o = { ...DefaultViewerOptions, ...options };

        const spec: PluginSpec = {
            actions: [...DefaultPluginSpec.actions],
            behaviors: [
                PluginSpec.Behavior(PluginBehaviors.Representation.HighlightLoci, { mark: false }),
                PluginSpec.Behavior(PluginBehaviors.Representation.DefaultLociLabelProvider),
                PluginSpec.Behavior(PluginBehaviors.Camera.FocusLoci),

                PluginSpec.Behavior(PluginBehaviors.CustomProps.StructureInfo),
                PluginSpec.Behavior(PluginBehaviors.CustomProps.Interactions),
                PluginSpec.Behavior(PluginBehaviors.CustomProps.SecondaryStructure),
            ],
            animations: [...DefaultPluginSpec.animations || []],
            customParamEditors: DefaultPluginSpec.customParamEditors,
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
                viewport: {
                    view: ViewportComponent
                }
            },
            config: [
                [PluginConfig.Viewport.ShowExpand, o.viewportShowExpand],
                [PluginConfig.Viewport.ShowControls, o.viewportShowControls],
                [PluginConfig.Viewport.ShowSettings, o.viewportShowSettings],
                [PluginConfig.Viewport.ShowSelectionMode, o.viewportShowSelectionMode],
                [PluginConfig.Viewport.ShowAnimation, o.viewportShowAnimation],
                [PluginConfig.State.DefaultServer, o.pluginStateServer],
                [PluginConfig.State.CurrentServer, o.pluginStateServer],
                [PluginConfig.VolumeStreaming.DefaultServer, o.volumeStreamingServer],
                [PluginConfig.Download.DefaultPdbProvider, o.pdbProvider],
                [PluginConfig.Download.DefaultEmdbProvider, o.emdbProvider]
            ]
        };

        const element = typeof elementOrId === 'string'
            ? document.getElementById(elementOrId)
            : elementOrId;
        if (!element) throw new Error(`Could not get element with id '${elementOrId}'`);
        this.plugin = createPlugin(element, spec);

        PluginCommands.Canvas3D.SetSettings(this.plugin, { settings: {
            renderer: {
                ...this.plugin.canvas3d!.props.renderer,
                backgroundColor: ColorNames.white,
            },
            camera: {
                ...this.plugin.canvas3d!.props.camera,
                helper: { axes: { name: 'off', params: {} } }
            }
        } });
    }

    setRemoteSnapshot(id: string) {
        const url = `${this.plugin.config.get(PluginConfig.State.CurrentServer)}/get/${id}`;
        return PluginCommands.State.Snapshots.Fetch(this.plugin, { url });
    }

    loadSnapshotFromUrl(url: string, type: PluginState.SnapshotType) {
        return PluginCommands.State.Snapshots.OpenUrl(this.plugin, { url, type });
    }

    loadStructureFromUrl(url: string, format: BuiltInTrajectoryFormat = 'mmcif', isBinary = false) {
        const params = DownloadStructure.createDefaultParams(this.plugin.state.data.root.obj!, this.plugin);
        return this.plugin.runTask(this.plugin.state.data.applyAction(DownloadStructure, {
            source: {
                name: 'url',
                params: {
                    url: Asset.Url(url),
                    format: format as any,
                    isBinary,
                    options: params.source.params.options,
                }
            }
        }));
    }

    async loadStructuresFromUrlsAndMerge(sources: { url: string, format: BuiltInTrajectoryFormat, isBinary?: boolean }[]) {
        const structures: { ref: string }[] = [];
        for (const { url, format, isBinary } of sources) {
            const data = await this.plugin.builders.data.download({ url, isBinary });
            const trajectory = await this.plugin.builders.structure.parseTrajectory(data, format);
            const model = await this.plugin.builders.structure.createModel(trajectory);
            const modelProperties = await this.plugin.builders.structure.insertModelProperties(model);
            const structure = await this.plugin.builders.structure.createStructure(modelProperties || model);
            const structureProperties = await this.plugin.builders.structure.insertStructureProperties(structure);

            structures.push({ ref: structureProperties?.ref || structure.ref });
        }
        const dependsOn = structures.map(({ ref }) => ref);
        const data = this.plugin.state.data.build().toRoot().apply(MergeStructures, { structures }, { dependsOn });
        const structure = await data.commit();
        const structureProperties = await this.plugin.builders.structure.insertStructureProperties(structure);
        await this.plugin.builders.structure.representation.applyPreset(structureProperties || structure, StructurePreset);
    }

    async loadStructureFromData(data: string | number[], format: BuiltInTrajectoryFormat, options?: { dataLabel?: string }) {
        const _data = await this.plugin.builders.data.rawData({ data, label: options?.dataLabel });
        const trajectory = await this.plugin.builders.structure.parseTrajectory(_data, format);
        await this.plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    }

    loadPdb(pdb: string) {
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
                    options: params.source.params.options,
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

    loadEmdb(emdb: string) {
        const provider = this.plugin.config.get(PluginConfig.Download.DefaultEmdbProvider)!;
        return this.plugin.runTask(this.plugin.state.data.applyAction(DownloadDensity, {
            source: {
                name: 'pdb-emd-ds' as const,
                params: {
                    provider: {
                        id: emdb,
                        server: provider,
                    },
                    detail: 3,
                }
            }
        }));
    }
}

type MergeStructures = typeof MergeStructures
const MergeStructures = PluginStateTransform.BuiltIn({
    name: 'merge-structures',
    display: { name: 'Merge Structures', description: 'Merge Structure' },
    from: PSO.Root,
    to: PSO.Molecule.Structure,
    params: {
        structures: PD.ObjectList({
            ref: PD.Text('')
        }, ({ ref }) => ref, { isHidden: true })
    }
})({
    apply({ params, dependencies }) {
        return Task.create('Merge Structures', async ctx => {
            if (params.structures.length === 0) return StateObject.Null;

            const first = dependencies![params.structures[0].ref].data as Structure;
            const builder = Structure.Builder({ masterModel: first.models[0] });
            for (const { ref } of params.structures) {
                const s = dependencies![ref].data as Structure;
                for (const unit of s.units) {
                    // TODO invariantId
                    builder.addUnit(unit.kind, unit.model, unit.conformation.operator, unit.elements, unit.traits);
                }
            }

            const structure = builder.getStructure();
            return new PSO.Molecule.Structure(structure, { label: 'Merged Structure' });
        });
    }
});

(window as any).DockingViewer = Viewer;