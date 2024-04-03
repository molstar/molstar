import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { PluginContext } from '../../../mol-plugin/context';
import { StateObjectRef } from '../../../mol-state';
import { GEOMETRIC_SEGMENTATION_NODE_TAG, MESH_SEGMENTATION_NODE_TAG, SEGMENTATION_NODE_TAG, VOLUME_NODE_TAG } from '../new-volumes-and-segmentations/entry-root';
import { ProjectGeometricSegmentationData, ProjectGeometricSegmentationDataParamsValues, ProjectLatticeSegmentationDataParamsValues, ProjectMeshData, ProjectMeshSegmentationDataParamsValues, ProjectSegmentationData, ProjectVolumeData, VolsegEntryFromFile, VolsegGlobalStateFromFile, VolsegStateFromEntry } from '../new-volumes-and-segmentations/transformers';
import { getSegmentLabelsFromDescriptions } from '../new-volumes-and-segmentations/volseg-api/utils';

export async function loadCVSXFromAnything(plugin: PluginContext, data: StateObjectRef<PluginStateObject.Data.Binary | PluginStateObject.Data.String>) {
    await plugin.build().to(data).apply(VolsegGlobalStateFromFile, {}, { state: { isGhost: true } }).commit();

    const entryNode = await plugin.build().to(data).apply(VolsegEntryFromFile).commit();
    await plugin.build().to(entryNode).apply(VolsegStateFromEntry, {}, { state: { isGhost: true } }).commit();


    if (entryNode.data) {
        const entryData = entryNode.data;
        let timeframeIndex = entryData.metadata!.value!.raw.grid.volumes.time_info.start;
        let channelIds = entryData.metadata!.value!.raw.grid.volumes.channel_ids;
        // const currentTimeframe = entryData.currentTimeframe.value;
        if (entryData.filesData!.query.channel_id) {
            channelIds = [entryData.filesData!.query.channel_id];
        }
        if (entryData.filesData!.query.time) {
            timeframeIndex = entryData.filesData!.query.time;
        }
        const hasVolumes = entryNode.data.metadata.value!.raw.grid.volumes.volume_sampling_info.spatial_downsampling_levels.length > 0;
        if (hasVolumes) {
            const group = await entryNode.data.volumeData.createVolumeGroup();
            const updatedChannelsData = [];
            const results = [];
            for (const channelId of channelIds) {
                const volumeParams = { timeframeIndex: timeframeIndex, channelId: channelId };
                const volumeNode = await plugin.build().to(group).apply(ProjectVolumeData, volumeParams, { tags: [VOLUME_NODE_TAG] }).commit();
                const result = await entryNode.data.volumeData.createVolumeRepresentation3D(volumeNode, volumeParams);
                results.push(result);
            }
            for (const result of results) {
                if (result) {
                    const isovalue = result.isovalue.kind === 'relative' ? result.isovalue.relativeValue : result.isovalue.absoluteValue;
                    updatedChannelsData.push(
                        {
                            channelId: result.channelId, volumeIsovalueKind: result.isovalue.kind, volumeIsovalueValue: isovalue, volumeType: result.volumeType, volumeOpacity: result.opacity,
                            label: result.label,
                            color: result.color
                        }
                    );
                }
            }
            await entryNode.data.updateStateNode({ channelsData: [...updatedChannelsData] });
        }

        const hasLattices = entryData.metadata!.value!.hasLatticeSegmentations();
        if (hasLattices) {
            let segmentationIds = hasLattices.segmentation_ids;
            if (entryData.filesData!.query.segmentation_id) {
                segmentationIds = [entryData.filesData!.query.segmentation_id];
            }
            // loop over lattices and create one for each
            const group = await entryNode.data.latticeSegmentationData.createSegmentationGroup();
            for (const segmentationId of segmentationIds) {
                // const segmentationId = hasLattices.segmentation_sampling_info
                // same, single channel single timeframe

                // const segmentationIds = hasLattices.segmentation_ids;
                // for (const segmentationId of segmentationIds) {
                const descriptionsForLattice = entryNode.data.metadata.value!.getAllDescriptionsForSegmentationAndTimeframe(
                    segmentationId,
                    'lattice',
                    timeframeIndex
                );
                const segmentLabels = getSegmentLabelsFromDescriptions(descriptionsForLattice);
                const segmentationParams: ProjectLatticeSegmentationDataParamsValues = {
                    timeframeIndex: timeframeIndex,
                    segmentationId: segmentationId,
                    segmentLabels: segmentLabels,
                    ownerId: entryNode.data.ref
                };
                const segmentationNode = await plugin.build().to(group).apply(ProjectSegmentationData, segmentationParams, { tags: [SEGMENTATION_NODE_TAG] }).commit();
                await entryNode.data.latticeSegmentationData.createSegmentationRepresentation3D(segmentationNode, segmentationParams);
            }

        }

        const hasGeometricSegmentation = entryData.metadata!.value!.hasGeometricSegmentations();
        if (hasGeometricSegmentation) {
            let segmentationIds = hasGeometricSegmentation.segmentation_ids;
            if (entryData.filesData!.query.segmentation_id) {
                segmentationIds = [entryData.filesData!.query.segmentation_id];
            }

            const group = await entryNode.data.geometricSegmentationData.createGeometricSegmentationGroup();
            for (const segmentationId of segmentationIds) {
                // const segmentationId: string = entryData.filesData!.query.args.segmentation_id;
                const geometricSegmentationParams: ProjectGeometricSegmentationDataParamsValues = {
                    segmentationId: segmentationId,
                    timeframeIndex: timeframeIndex
                };
                const geometricSegmentationNode = await plugin.build().to(group).apply(ProjectGeometricSegmentationData, geometricSegmentationParams, { tags: [GEOMETRIC_SEGMENTATION_NODE_TAG] }).commit();
                await entryNode.data.geometricSegmentationData.createGeometricSegmentationRepresentation3D(geometricSegmentationNode, geometricSegmentationParams);
            // }
            }
        }

        const hasMeshes = entryNode.data.metadata.value!.hasMeshSegmentations();
        if (hasMeshes) {
            let segmentationIds = hasMeshes.segmentation_ids;
            if (entryData.filesData!.query.segmentation_id) {
                segmentationIds = [entryData.filesData!.query.segmentation_id];
            }
            const group = await entryNode.data.meshSegmentationData.createMeshGroup();
            for (const segmentationId of segmentationIds) {
                const meshSegmentParams = entryData.meshSegmentationData.getMeshSegmentParams(segmentationId, timeframeIndex);
                const meshParams: ProjectMeshSegmentationDataParamsValues = {
                    meshSegmentParams: meshSegmentParams,
                    segmentationId: segmentationId,
                    timeframeIndex: timeframeIndex
                };
                const meshNode = await plugin.build().to(group).apply(ProjectMeshData, meshParams, { tags: [MESH_SEGMENTATION_NODE_TAG] }).commit();
                await entryNode.data.meshSegmentationData.createMeshRepresentation3D(meshNode, meshParams);
            }
        }


    };
    return entryNode;
}