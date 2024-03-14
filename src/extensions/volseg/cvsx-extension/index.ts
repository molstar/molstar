import { PluginStateObject } from "../../../mol-plugin-state/objects";
import { PluginContext } from "../../../mol-plugin/context";
import { StateObjectSelector } from "../../../mol-state";
import { GEOMETRIC_SEGMENTATION_NODE_TAG, MESH_SEGMENTATION_NODE_TAG, SEGMENTATION_NODE_TAG, VOLUME_NODE_TAG } from "../new-volumes-and-segmentations/entry-root";
import { ProjectGeometricSegmentationData, ProjectGeometricSegmentationDataParamsValues, ProjectLatticeSegmentationDataParamsValues, ProjectMeshData, ProjectMeshSegmentationDataParamsValues, ProjectSegmentationData, ProjectVolumeData, VolsegEntryFromFile, VolsegGlobalStateFromFile, VolsegGlobalStateFromRoot, VolsegStateFromEntry } from "../new-volumes-and-segmentations/transformers";
import { getSegmentLabelsFromDescriptions } from "../new-volumes-and-segmentations/volseg-api/utils";

export async function loadCVSXFromAnything(plugin: PluginContext, data: StateObjectSelector<PluginStateObject.Data.Binary | PluginStateObject.Data.String>) {
    await plugin.build().to(data).apply(VolsegGlobalStateFromFile, {}, { state: { isGhost: true } }).commit();

    const entryNode = await plugin.build().to(data).apply(VolsegEntryFromFile).commit();
    await plugin.build().to(entryNode).apply(VolsegStateFromEntry, {}, { state: { isGhost: true } }).commit();
    debugger;

    if (entryNode.data) {
        const entryData = entryNode.data;
        const timeframeIndex = entryData.filesData!.query.args.time;
        const channelId = entryData.filesData!.query.args.channel_id;
        // const currentTimeframe = entryData.currentTimeframe.value;
        const hasVolumes = entryNode.data.metadata.value!.raw.grid.volumes.volume_sampling_info.spatial_downsampling_levels.length > 0;
        if (hasVolumes) {
            const group = await entryNode.data.volumeData.createVolumeGroup();
            const updatedChannelsData = [];
            const results = [];

            // single channel, single timeframe, get from entryData
            // const channelIds = entryNode.data.metadata.value!.raw.grid.volumes.channel_ids;
            // for (const channelId of channelIds) {
            const volumeParams = { timeframeIndex: timeframeIndex, channelId: channelId };
            const volumeNode = await plugin.build().to(group).apply(ProjectVolumeData, volumeParams, { tags: [VOLUME_NODE_TAG] }).commit();
            const result = await entryNode.data.volumeData.createVolumeRepresentation3D(volumeNode, volumeParams);
            results.push(result);
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

        const hasLattices = entryNode.data.metadata.value!.raw.grid.segmentation_lattices;
        if (hasLattices && hasLattices.segmentation_ids.length > 0) {
            const segmentationId: string = entryData.filesData!.query.args.segmentation_id;
            const group = await entryNode.data.latticeSegmentationData.createSegmentationGroup();
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

        const hasGeometricSegmentation = entryData.metadata.value!.raw.grid.geometric_segmentation;
        if (hasGeometricSegmentation && hasGeometricSegmentation.segmentation_ids.length > 0) {
            const group = await entryNode.data.geometricSegmentationData.createGeometricSegmentationGroup();
            // const timeInfo = entryData.metadata.value!.raw.grid.geometric_segmentation!.time_info;
            // single segmentation id
            // for (const segmentationId of hasGeometricSegmentation.segmentation_ids) {
            // const timeframeIndex = 0;
            const segmentationId: string = entryData.filesData!.query.args.segmentation_id;
            const geometricSegmentationParams: ProjectGeometricSegmentationDataParamsValues = {
                segmentationId: segmentationId,
                timeframeIndex: timeframeIndex
            };
            const geometricSegmentationNode = await plugin.build().to(group).apply(ProjectGeometricSegmentationData, geometricSegmentationParams, { tags: [GEOMETRIC_SEGMENTATION_NODE_TAG] }).commit();
            await entryNode.data.geometricSegmentationData.createGeometricSegmentationRepresentation3D(geometricSegmentationNode, geometricSegmentationParams);
            // }
        }

        const hasMeshes = entryNode.data.metadata.value!.raw.grid.segmentation_meshes;
        if (hasMeshes && hasMeshes.segmentation_ids.length > 0) {
            // meshes should be rendered as segmentation sets similar to lattices
            const group = await entryNode.data.meshSegmentationData.createMeshGroup();
            // const segmentationIds = hasMeshes.segmentation_ids;
            // for (const segmentationId of segmentationIds) {
            // const timeframeIndex = timeframeIndex;
            const segmentationId: string = entryData.filesData!.query.args.segmentation_id;
            const meshSegmentParams = entryData.meshSegmentationData.getMeshSegmentParams(segmentationId, timeframeIndex);
            const meshParams: ProjectMeshSegmentationDataParamsValues = {
                meshSegmentParams: meshSegmentParams,
                segmentationId: segmentationId,
                timeframeIndex: timeframeIndex
            };
            const meshNode = await plugin.build().to(group).apply(ProjectMeshData, meshParams, { tags: [MESH_SEGMENTATION_NODE_TAG] }).commit();
            await entryNode.data.meshSegmentationData.createMeshRepresentation3D(meshNode, meshParams);
        }


    };
    return entryNode;
}