/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { PluginCommands } from '../../mol-plugin/commands';
import { StateObjectSelector } from '../../mol-state';
import { VolsegEntryData } from './entry-root';
import { CreateShapePrimitiveProvider, VolsegShapePrimitivesData } from './shape_primitives';
import { ProjectGeometricSegmentationDataParamsValues, VolsegGeometricSegmentation } from './transformers';
import { GeometricSegmentationData } from './volseg-api/data';


const GEOMETRIC_SEGMENTATION_GROUP_TAG = 'geometric-segmentation-group';

export class VolsegGeometricSegmentationData {
    private entryData: VolsegEntryData;

    constructor(rootData: VolsegEntryData) {
        this.entryData = rootData;
    }

    async createGeometricSegmentationGroup() {
        let group = this.entryData.findNodesByTags(GEOMETRIC_SEGMENTATION_GROUP_TAG)[0]?.transform.ref;
        if (!group) {
            const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup,
                { label: 'Segmentation', description: 'Geometric segmentation' }, { tags: [GEOMETRIC_SEGMENTATION_GROUP_TAG], state: { isCollapsed: true } }).commit();
            group = newGroupNode.ref;
        }
        return group;
    }

    async createGeometricSegmentationRepresentation3D(gsNode: StateObjectSelector<VolsegGeometricSegmentation>, params: ProjectGeometricSegmentationDataParamsValues) {
        const gsData: VolsegShapePrimitivesData = gsNode.cell!.obj!.data;
        // const update = this.entryData.plugin.build().to(gsNode);
        const { timeframeIndex, segmentationId } = params;
        const descriptions = this.entryData.metadata.value!.getAllDescriptionsForSegmentationAndTimeframe(segmentationId, 'primitive', timeframeIndex);
        const segmentAnnotations = this.entryData.metadata.value!.getAllSegmentAnotationsForSegmentationAndTimeframe(segmentationId, 'primitive', timeframeIndex);
        for (const primitiveData of gsData.shapePrimitiveData.shape_primitive_list) {
            await this.entryData.newUpdate().to(gsNode.ref)
            // TODO: can provide a single description and a single segment annotation
                .apply(CreateShapePrimitiveProvider, { segmentId: primitiveData.id, descriptions: descriptions, segmentAnnotations: segmentAnnotations, segmentationId: segmentationId })
                // TODO: shape representation 3d could have no alpha
                .apply(StateTransforms.Representation.ShapeRepresentation3D, { alpha: 0.5 }, { tags: ['geometric-segmentation-visual', segmentationId, `segment-${primitiveData.id}`] })
                .commit();
        }
        // await update.commit();
    }

    // async loadGeometricSegmentation(timeframeIndex: number) {
    //     const hasGeometricSegmentation = this.entryData.metadata.value!.raw.grid.geometric_segmentation;
    //     if (hasGeometricSegmentation && hasGeometricSegmentation.segmentation_ids.length > 0) {
    //         let group = this.entryData.findNodesByTags(GEOMETRIC_SEGMENTATION_GROUP_TAG)[0]?.transform.ref;
    //         if (!group) {
    //             const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup,
    //                 { label: 'Segmentation', description: 'Geometric segmentation' }, { tags: [GEOMETRIC_SEGMENTATION_GROUP_TAG], state: { isCollapsed: true } }).commit();
    //             group = newGroupNode.ref;
    //         }
    //         // const timeInfo = this.entryData.metadata.value!.raw.grid.geometric_segmentation!.time_info;
    //         for (const segmentationId of hasGeometricSegmentation.segmentation_ids) {
    //             const url = this.entryData.api.geometricSegmentationUrl(this.entryData.source, this.entryData.entryId, segmentationId);

    //             const primitivesData = await this.entryData._resolveStringUrl(url);

    //             const parsedData: GeometricSegmentationData = JSON.parse(primitivesData);
    //             console.log('parsedData', parsedData);
    //             // const t = timeInfo[segmentationId];
    //             // for (let timeframeIndex = t.start; timeframeIndex <= t.end; timeframeIndex++) {
    //             const timeframeData = parsedData.primitives[timeframeIndex];
    //             const descriptions = this.entryData.metadata.value!.getAllDescriptionsForSegmentationAndTimeframe(segmentationId, 'primitive', timeframeIndex);
    //             const segmentAnnotations = this.entryData.metadata.value!.getAllSegmentAnotationsForSegmentationAndTimeframe(segmentationId, 'primitive', timeframeIndex);
    //             for (const shapePrimitiveData of timeframeData.shape_primitive_list) {
    //                 const shapePrimitiveNode = await this.entryData.newUpdate().to(group)
    //                 // TODO: can provide a single description and a single segment annotation
    //                     .apply(CreateShapePrimitiveProvider, { data: shapePrimitiveData, descriptions: descriptions, segmentAnnotations: segmentAnnotations, segmentationId: segmentationId })
    //                     // TODO: shape representation 3d could have no alpha
    //                     .apply(StateTransforms.Representation.ShapeRepresentation3D, { alpha: 0.5 }, { tags: ['geometric-segmentation-visual', segmentationId, `segment-${shapePrimitiveData.id}`] })
    //                     .commit();
    //             }
    //         }
    //     }
    // }
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

    async selectSegment(segment?: number, segmentationId?: string) {
        if (segment === undefined || segment < 0 || segmentationId === undefined) return;
        const visuals = this.entryData.findNodesByTags('geometric-segmentation-visual', `segment-${segment}`, segmentationId);
        const reprNode: PluginStateObject.Shape.Representation3D | undefined = visuals[0]?.obj;
        debugger;
        if (!reprNode) return;
        const loci = reprNode.data.repr.getAllLoci()[0];
        debugger;
        if (!loci) return;
        this.entryData.plugin.managers.interactivity.lociSelects.select({ loci: loci, repr: reprNode.data.repr }, false);
    }
}