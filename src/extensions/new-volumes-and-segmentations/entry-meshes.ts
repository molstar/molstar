/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginStateObject } from '../../mol-plugin-state/objects';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { ShapeRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { PluginCommands } from '../../mol-plugin/commands';
import { Color } from '../../mol-util/color';
import { ColorNames } from '../../mol-util/color/names';
import { Box3D } from '../../mol-math/geometry';

import { BACKGROUND_SEGMENT_VOLUME_THRESHOLD } from '../meshes/mesh-streaming/behavior';
import { BACKGROUND_OPACITY, CreateMeshlistStateObject, FOREROUND_OPACITY, MeshShapeTransformer, MeshlistData, VolsegMeshSegmentation, createMeshFromUrl, meshSegmentParamsValues } from '../meshes/mesh-extension';
import { Segment } from './volseg-api/data';
import { VolsegEntryData } from './entry-root';


import { StateObjectSelector } from '../../mol-state';

export const DEFAULT_MESH_DETAIL: number | null = 5; // null means worst

export class VolsegMeshSegmentationData {
    private entryData: VolsegEntryData;

    constructor(rootData: VolsegEntryData) {
        this.entryData = rootData;
    }

    async loadSegmentation() {
        const hasMeshes = this.entryData.metadata.meshSegmentIds.length > 0;
        if (hasMeshes) {
            await this.showSegments(this.entryData.metadata.allSegmentIds);
        }
    }

    updateOpacity(opacity: number) {
        const visuals = this.entryData.findNodesByTags('mesh-segment-visual');
        const update = this.entryData.newUpdate();
        for (const visual of visuals) {
            update.to(visual).update(ShapeRepresentation3D, p => { (p as any).alpha = opacity; });
        }
        return update.commit();
    }

    async createMeshRepresentation3D(meshNode: StateObjectSelector<VolsegMeshSegmentation>) {
        const meshData = meshNode.data!.meshData;
        const totalVolume = this.entryData.metadata.gridTotalVolume;
        for (const meshDataItem of meshData) {
            const update = this.entryData.plugin.build().to(meshNode);
            const meshListStateObj = await update.to(meshNode)
                .apply(CreateMeshlistStateObject,
                    { segmentId: meshDataItem.meshSegmentParams.id }
                )
                .commit();

            const transparentIfBboxAbove = BACKGROUND_SEGMENT_VOLUME_THRESHOLD * totalVolume;

            let transparent = false;
            if (transparentIfBboxAbove !== undefined && meshListStateObj.data) {
                const bbox = MeshlistData.bbox(meshListStateObj.data) || Box3D.zero();
                transparent = Box3D.volume(bbox) > transparentIfBboxAbove;
            }

            await this.entryData.plugin.build().to(meshListStateObj)
                .apply(MeshShapeTransformer, { color: meshDataItem.meshSegmentParams.color },)
                .apply(ShapeRepresentation3D,
                    { alpha: transparent ? BACKGROUND_OPACITY : FOREROUND_OPACITY },
                    { tags: ['mesh-segment-visual', `segment-${meshDataItem.meshSegmentParams.id}`] }
                )
                .commit();
        }
        this.entryData.actionShowSegments(this.entryData.metadata.meshSegmentIds);
    }

    getMeshSegmentParams() {
        const params: meshSegmentParamsValues[] = [];
        const segmentsToCreate = this.entryData.metadata.meshSegmentIds;
        for (const seg of segmentsToCreate) {
            const segment = this.entryData.metadata.getSegment(seg);
            if (!segment) continue;
            const detail = this.entryData.metadata.getSufficientMeshDetail(seg, DEFAULT_MESH_DETAIL);
            const color = segment.color.length >= 3 ? Color.fromNormalizedArray(segment.color, 0) : ColorNames.gray;
            // const url = this.entryData.api.meshUrl_Bcif(this.entryData.source, this.entryData.entryId, seg, detail, timeframeIndex, channelId);
            const label = segment.biological_annotation.name ?? `Segment ${seg}`;
            const segmentParams: meshSegmentParamsValues = {
                id: seg,
                label: label,
                // url: url,
                color: color,
                detail: detail
            };
            params.push(segmentParams)
        }

        return params;
    }

    async createMeshGroup() {
        let group = this.entryData.findNodesByTags('mesh-segmentation-group')[0]?.transform.ref;
        if (!group) {
            const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup, { label: 'Segmentation', description: 'Mesh' }, { tags: ['mesh-segmentation-group'], state: { isCollapsed: true } }).commit();
            group = newGroupNode.ref;
        }
        return group;
    }

    async highlightSegment(segment: Segment) {
        const visuals = this.entryData.findNodesByTags('mesh-segment-visual', `segment-${segment.id}`);
        for (const visual of visuals) {
            await PluginCommands.Interactivity.Object.Highlight(this.entryData.plugin, { state: this.entryData.plugin.state.data, ref: visual.transform.ref });
        }
    }

    async selectSegment(segment?: number) {
        if (segment === undefined || segment < 0) return;
        const visuals = this.entryData.findNodesByTags('mesh-segment-visual', `segment-${segment}`);
        const reprNode: PluginStateObject.Shape.Representation3D | undefined = visuals[0]?.obj;
        if (!reprNode) return;
        const loci = reprNode.data.repr.getAllLoci()[0];
        if (!loci) return;
        this.entryData.plugin.managers.interactivity.lociSelects.select({ loci: loci, repr: reprNode.data.repr }, false);
    }

    /** Make visible the specified set of mesh segments */
    async showSegments(segments: number[]) {
        const segmentsToShow = new Set(segments);

        const visuals = this.entryData.findNodesByTags('mesh-segment-visual');
        for (const visual of visuals) {
            const theTag = visual.obj?.tags?.find(tag => tag.startsWith('segment-'));
            if (!theTag) continue;
            const id = parseInt(theTag.split('-')[1]);
            const visibility = segmentsToShow.has(id);
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.transform.ref, !visibility); // true means hide, ¯\_(ツ)_/¯
            segmentsToShow.delete(id);
        }

        const segmentsToCreate = this.entryData.metadata.meshSegmentIds.filter(seg => segmentsToShow.has(seg));
        if (segmentsToCreate.length === 0) return;

        let group = this.entryData.findNodesByTags('mesh-segmentation-group')[0]?.transform.ref;
        if (!group) {
            const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup, { label: 'Segmentation', description: 'Mesh' }, { tags: ['mesh-segmentation-group'], state: { isCollapsed: true } }).commit();
            group = newGroupNode.ref;
        }
        const totalVolume = this.entryData.metadata.gridTotalVolume;

        const awaiting = [];
        for (const seg of segmentsToCreate) {
            const segment = this.entryData.metadata.getSegment(seg);
            if (!segment) continue;
            const detail = this.entryData.metadata.getSufficientMeshDetail(seg, DEFAULT_MESH_DETAIL);
            const color = segment.color.length >= 3 ? Color.fromNormalizedArray(segment.color, 0) : ColorNames.gray;
            const url = this.entryData.api.meshUrl_Bcif(this.entryData.source, this.entryData.entryId, seg, detail, this.entryData.currentTimeframe.value, 0);
            const label = segment.biological_annotation.name ?? `Segment ${seg}`;
            const meshPromise = createMeshFromUrl(this.entryData.plugin, url, seg, detail, true, color, group,
                BACKGROUND_SEGMENT_VOLUME_THRESHOLD * totalVolume, `<b>${label}</b>`, this.entryData.ref);
            awaiting.push(meshPromise);
        }
        for (const promise of awaiting) await promise;
    }
}
