/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { StateTransforms } from '../../mol-plugin-state/transforms';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { PluginCommands } from '../../mol-plugin/commands';
import { VolsegEntryData } from './entry-root';
import { CreateShapePrimitiveProvider } from './shape_primitives';
import { GeometricSegmentationData } from './volseg-api/data';


const GEOMETRIC_SEGMENTATION_GROUP_TAG = 'geometric-segmentation-group';

export class VolsegGeometricSegmentationData {
    private entryData: VolsegEntryData;

    constructor(rootData: VolsegEntryData) {
        this.entryData = rootData;
    }

    async loadGeometricSegmentation(timeframeIndex: number) {
        const hasGeometricSegmentation = this.entryData.metadata.raw.grid.geometric_segmentation;
        if (hasGeometricSegmentation && hasGeometricSegmentation.segmentation_ids.length > 0) {
            let group = this.entryData.findNodesByTags(GEOMETRIC_SEGMENTATION_GROUP_TAG)[0]?.transform.ref;
            if (!group) {
                const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup,
                    { label: 'Segmentation', description: 'Geometric segmentation' }, { tags: [GEOMETRIC_SEGMENTATION_GROUP_TAG], state: { isCollapsed: true } }).commit();
                group = newGroupNode.ref;
            }
            // const timeInfo = this.entryData.metadata.raw.grid.geometric_segmentation!.time_info;
            for (const segmentationId of hasGeometricSegmentation.segmentation_ids) {
                const url = this.entryData.api.geometricSegmentationUrl(this.entryData.source, this.entryData.entryId, segmentationId);

                const primitivesData = await this.entryData._resolveStringUrl(url);

                const parsedData: GeometricSegmentationData = JSON.parse(primitivesData);
                console.log('parsedData', parsedData);
                // const t = timeInfo[segmentationId];
                // for (let timeframeIndex = t.start; timeframeIndex <= t.end; timeframeIndex++) {
                const timeframeData = parsedData.primitives[timeframeIndex];
                const descriptions = this.entryData.metadata.getAllDescriptionsForSegmentationAndTimeframe(segmentationId, 'primitive', timeframeIndex);
                const segmentAnnotations = this.entryData.metadata.getAllSegmentAnotationsForSegmentationAndTimeframe(segmentationId, 'primitive', timeframeIndex);
                for (const shapePrimitiveData of timeframeData.shape_primitive_list) {
                    const shapePrimitiveNode = await this.entryData.newUpdate().to(group)
                    // TODO: can provide a single description and a single segment annotation
                        .apply(CreateShapePrimitiveProvider, { data: shapePrimitiveData, descriptions: descriptions, segmentAnnotations: segmentAnnotations })
                        // TODO: shape representation 3d could have no alpha
                        .apply(StateTransforms.Representation.ShapeRepresentation3D, { alpha: 0.5 }, { tags: ['geometric-segmentation-visual', segmentationId, `segment-${shapePrimitiveData.id}`] })
                        .commit();
                }
            }
        }
    }
    // From meshes, here probably similar
    async showSegments(segmentIds: number[], segmentationId: string) {
        const segmentsToShow = new Set(segmentIds);

        // This will select all segments of that segmentation
        const visuals = this.entryData.findNodesByTags('geometric-segmentation-visual', segmentationId);
        for (const visual of visuals) {
            const theTag = visual.obj?.tags?.find(tag => tag.startsWith('segment-'));
            if (!theTag) continue;
            const id = parseInt(theTag.split('-')[1]);
            const visibility = segmentsToShow.has(id);
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.transform.ref, !visibility); // true means hide, ¯\_(ツ)_/¯
            segmentsToShow.delete(id);
        }
    }

    async highlightSegment(segmentId: number, segmentationId: string) {
        const visuals = this.entryData.findNodesByTags('geometric-segmentation-visual', `segment-${segmentId}`, segmentationId);
        for (const visual of visuals) {
            await PluginCommands.Interactivity.Object.Highlight(this.entryData.plugin, { state: this.entryData.plugin.state.data, ref: visual.transform.ref });
        }
    }
}