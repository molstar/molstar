/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginStateObject as SO } from '../../mol-plugin-state/objects';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginConfigItem } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { StateAction } from '../../mol-state';
import { Task } from '../../mol-task';
import { DEFAULT_VOLSEG_SERVER, VolumeApiV2 } from './volseg-api/api';

import { GEOMETRIC_SEGMENTATION_NODE_TAG, SEGMENTATION_NODE_TAG, VOLUME_NODE_TAG, VolsegEntryData, VolsegEntryParamValues, createLoadVolsegParams } from './entry-root';
import { VolsegGlobalState } from './global-state';
import { createEntryId } from './helpers';
import { ProjectGeometricSegmentationData, ProjectGeometricSegmentationDataParamsValues, ProjectMeshData, ProjectMeshSegmentationDataParamsValues, ProjectSegmentationData, ProjectSegmentationDataParamsValues, ProjectVolumeData, VolsegEntryFromRoot, VolsegGlobalStateFromRoot, VolsegStateFromEntry } from './transformers';
import { VolsegUI } from './ui';
import { createSegmentKey, getSegmentLabelsFromDescriptions } from './volseg-api/utils';
import { useBehavior } from '../../mol-plugin-ui/hooks/use-behavior';

// TODO: temp change, put there 'localhost'
const DEBUGGING = typeof window !== 'undefined' ? window?.location?.hostname === 'localhost' || '127.0.0.1' : false;

export const NewVolsegVolumeServerConfig = {
    // DefaultServer: new PluginConfigItem('volseg-volume-server', DEFAULT_VOLUME_SERVER_V2),
    DefaultServer: new PluginConfigItem('volseg-volume-server', DEBUGGING ? 'http://localhost:9000/v2' : DEFAULT_VOLSEG_SERVER),
};


export const NewVolseg = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'new-volseg',
    category: 'misc',
    display: {
        name: 'New Volseg',
        description: 'New Volseg'
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        register() {
            this.ctx.state.data.actions.add(LoadVolseg);
            this.ctx.customStructureControls.set('new-volseg', VolsegUI as any);
            this.initializeEntryLists(); // do not await

            const entries = new Map<string, VolsegEntryData>();
            this.subscribeObservable(this.ctx.state.data.events.cell.created, o => {
                if (o.cell.obj instanceof VolsegEntryData) entries.set(o.ref, o.cell.obj);
            });

            this.subscribeObservable(this.ctx.state.data.events.cell.removed, o => {
                if (entries.has(o.ref)) {
                    entries.get(o.ref)!.dispose();
                    entries.delete(o.ref);
                }
            });
        }
        unregister() {
            this.ctx.state.data.actions.remove(LoadVolseg);
            this.ctx.customStructureControls.delete('new-volseg');
        }
        private async initializeEntryLists() {
            const apiUrl = this.ctx.config.get(NewVolsegVolumeServerConfig.DefaultServer) ?? DEFAULT_VOLSEG_SERVER;
            const api = new VolumeApiV2(apiUrl);
            const entryLists = await api.getEntryList(10 ** 6);
            Object.values(entryLists).forEach(l => l.sort());
            (this.ctx.customState as any).volsegAvailableEntries = entryLists;
        }
    }
});


export const LoadVolseg = StateAction.build({
    display: { name: 'Load New Volume & Segmentation' },
    from: SO.Root,
    params: (a, plugin: PluginContext) => {
        const res = createLoadVolsegParams(plugin, (plugin.customState as any).volsegAvailableEntries);
        return res;
    },
})(({ params, state }, ctx: PluginContext) => Task.create('Loading Volume & Segmentation', taskCtx => {
    return state.transaction(async () => {
        const entryParams = VolsegEntryParamValues.fromLoadVolsegParamValues(params);
        if (entryParams.entryId.trim().length === 0) {
            alert('Must specify Entry Id!');
            throw new Error('Specify Entry Id');
        }
        if (!entryParams.entryId.includes('-')) {
            // add source prefix if the user omitted it (e.g. 1832 -> emd-1832)
            entryParams.entryId = createEntryId(entryParams.source, entryParams.entryId);
        }
        ctx.behaviors.layout.leftPanelTabName.next('data');

        const globalStateNode = ctx.state.data.selectQ(q => q.ofType(VolsegGlobalState))[0];
        if (!globalStateNode) {
            await state.build().toRoot().apply(VolsegGlobalStateFromRoot, {}, { state: { isGhost: !DEBUGGING } }).commit();
        }

        const entryNode = await state.build().toRoot().apply(VolsegEntryFromRoot, entryParams).commit();
        await state.build().to(entryNode).apply(VolsegStateFromEntry, {}, { state: { isGhost: !DEBUGGING } }).commit();
        if (entryNode.data) {
            const entryData = entryNode.data;
            // const currentTimeframe = entryData.currentTimeframe.value;
            const hasVolumes = entryNode.data.metadata.raw.grid.volumes.volume_sampling_info.spatial_downsampling_levels.length > 0;
            if (hasVolumes) {
                const group = await entryNode.data.volumeData.createVolumeGroup();
                const updatedChannelsData = [];
                const results = [];
                const channelIds = entryNode.data.metadata.raw.grid.volumes.channel_ids;
                for (const channelId of channelIds) {
                    const volumeParams = { timeframeIndex: 0, channelId: channelId };
                    const volumeNode = await state.build().to(group).apply(ProjectVolumeData, volumeParams, { tags: [VOLUME_NODE_TAG] }).commit();
                    const result = await entryNode.data.volumeData.createVolumeRepresentation3D(volumeNode, volumeParams);
                    results.push(result);
                }
                for (const result of results) {
                    if (result) {
                        const isovalue = result.isovalue.kind === 'relative' ? result.isovalue.relativeValue : result.isovalue.absoluteValue;
                        updatedChannelsData.push(
                            { channelId: result.channelId, volumeIsovalueKind: result.isovalue.kind, volumeIsovalueValue: isovalue, volumeType: result.volumeType, volumeOpacity: result.opacity,
                                label: result.label,
                                color: result.color
                            }
                        );
                    }
                }
                await entryNode.data.updateStateNode({ channelsData: [...updatedChannelsData] });
            }

            const hasLattices = entryNode.data.metadata.raw.grid.segmentation_lattices;
            if (hasLattices && hasLattices.segmentation_ids.length > 0) {
                const group = await entryNode.data.latticeSegmentationData.createSegmentationGroup();
                const segmentationIds = hasLattices.segmentation_ids;
                for (const segmentationId of segmentationIds) {
                    const descriptionsForLattice = entryNode.data.metadata.getAllDescriptionsForSegmentationAndTimeframe(
                        segmentationId,
                        'lattice',
                        0
                    );
                    const segmentLabels = getSegmentLabelsFromDescriptions(descriptionsForLattice);
                    const segmentationParams: ProjectSegmentationDataParamsValues = {
                        timeframeIndex: 0,
                        segmentationId: segmentationId,
                        segmentLabels: segmentLabels,
                        ownerId: entryNode.data.ref
                    };
                    const segmentationNode = await state.build().to(group).apply(ProjectSegmentationData, segmentationParams, { tags: [SEGMENTATION_NODE_TAG] }).commit();
                    await entryNode.data.latticeSegmentationData.createSegmentationRepresentation3D(segmentationNode, segmentationParams);
                }
            };

            const hasMeshes = entryNode.data.metadata.raw.grid.segmentation_meshes;
            if (hasMeshes && hasMeshes.segmentation_ids.length > 0) {
                // meshes should be rendered as segmentation sets similar to lattices
                const group = await entryNode.data.meshSegmentationData.createMeshGroup();
                const segmentationIds = hasMeshes.segmentation_ids;
                for (const segmentationId of segmentationIds) {
                    const timeframeIndex = 0;
                    const meshSegmentParams = entryData.meshSegmentationData.getMeshSegmentParams(segmentationId, timeframeIndex);
                    const meshParams: ProjectMeshSegmentationDataParamsValues = {
                        meshSegmentParams: meshSegmentParams,
                        segmentationId: segmentationId,
                        timeframeIndex: timeframeIndex
                    };
                    const meshNode = await state.build().to(group).apply(ProjectMeshData, meshParams).commit();
                    await entryNode.data.meshSegmentationData.createMeshRepresentation3D(meshNode, meshParams);
                }


            }
            // for now for a single timeframe;
            // await entryData.geometricSegmentationData.loadGeometricSegmentation(0);
            const hasGeometricSegmentation = entryData.metadata.raw.grid.geometric_segmentation;
            if (hasGeometricSegmentation && hasGeometricSegmentation.segmentation_ids.length > 0) {
                const group = await entryNode.data.geometricSegmentationData.createGeometricSegmentationGroup();
                // const timeInfo = this.entryData.metadata.raw.grid.geometric_segmentation!.time_info;
                for (const segmentationId of hasGeometricSegmentation.segmentation_ids) {
                    const timeframeIndex = 0;
                    const geometricSegmentationParams: ProjectGeometricSegmentationDataParamsValues = {
                        segmentationId: segmentationId,
                        timeframeIndex: timeframeIndex
                    };
                    const geometricSegmentationNode = await state.build().to(group).apply(ProjectGeometricSegmentationData, geometricSegmentationParams, { tags: [GEOMETRIC_SEGMENTATION_NODE_TAG] }).commit();
                    await entryNode.data.geometricSegmentationData.createGeometricSegmentationRepresentation3D(geometricSegmentationNode, geometricSegmentationParams);
                }
            }
            const allAnnotationsForTimeframe = entryData.metadata.getAllAnnotationsForTimeframe(0);
            const allSegmentKeysForTimeframe = allAnnotationsForTimeframe.map(a => {
                return createSegmentKey(a.segment_id, a.segmentation_id, a.segment_kind);
            }
            );
            await entryData.actionShowSegments(allSegmentKeysForTimeframe);
        }
    }).runInContext(taskCtx);
}));