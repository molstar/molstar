/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Volume } from '../../mol-model/volume';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Download, ParseCif } from '../../mol-plugin-state/transforms/data';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { VolumeFromSegmentationCif } from '../../mol-plugin-state/transforms/volume';
import { PluginCommands } from '../../mol-plugin/commands';
import { Color } from '../../mol-util/color';

import { BOX, VolsegEntryData, MAX_VOXELS } from './entry-root';
import { VolumeVisualParams } from './entry-volume';
import { VolsegGlobalStateData } from './global-state';
import { StateObjectSelector } from '../../mol-state';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { ProjectSegmentationDataParamsValues } from './transformers';
import { createSegmentKey, parseSegmentKey } from './volseg-api/utils';


const GROUP_TAG = 'lattice-segmentation-group';
export const SEGMENT_VISUAL_TAG = 'lattice-segment-visual';

const DEFAULT_SEGMENT_COLOR = Color.fromNormalizedRgb(0.8, 0.8, 0.8);


export class VolsegLatticeSegmentationData {
    private entryData: VolsegEntryData;

    constructor(rootData: VolsegEntryData) {
        this.entryData = rootData;
    }

    async loadSegmentation() {
        // const hasLattices = this.entryData.metadata.raw.grid.segmentation_lattices.segmentation_ids.length > 0;
        // if (hasLattices) {
        //     const url = this.entryData.api.latticeUrl(this.entryData.source, this.entryData.entryId, 0, BOX, MAX_VOXELS);
        //     let group = this.entryData.findNodesByTags(GROUP_TAG)[0]?.transform.ref;
        //     if (!group) {
        //         const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup,
        //             { label: 'Segmentation', description: 'Lattice' }, { tags: [GROUP_TAG], state: { isCollapsed: true } }).commit();
        //         group = newGroupNode.ref;
        //     }
        //     const segmentLabels = this.entryData.metadata.allSegments.map(seg => ({ id: seg.id, label: seg.biological_annotation.name ? `<b>${seg.biological_annotation.name}</b>` : '' }));
        //     const volumeNode = await this.entryData.newUpdate().to(group)
        //         .apply(Download, { url, isBinary: true, label: `Segmentation Data: ${url}` })
        //         .apply(ParseCif)
        //         .apply(VolumeFromSegmentationCif, { blockHeader: 'SEGMENTATION_DATA', segmentLabels: segmentLabels, ownerId: this.entryData.ref })
        //         .commit();
        //     const volumeData = volumeNode.data as Volume;
        //     const segmentation = Volume.Segmentation.get(volumeData);
        //     const segmentIds: number[] = Array.from(segmentation?.segments.keys() ?? []);
        //     await this.entryData.newUpdate().to(volumeNode)
        //         .apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.entryData.plugin, volumeData, {
        //             type: 'segment',
        //             typeParams: { tryUseGpu: VolsegGlobalStateData.getGlobalState(this.entryData.plugin)?.tryUseGpu },
        //             color: 'volume-segment',
        //             colorParams: { palette: this.createPalette(segmentIds) },
        //         }), { tags: [SEGMENT_VISUAL_TAG] }).commit();
        // }
    }

    async createSegmentationGroup() {
        let group = this.entryData.findNodesByTags(GROUP_TAG)[0]?.transform.ref;
        if (!group) {
            const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup,
                { label: 'Segmentation', description: 'Lattice' }, { tags: [GROUP_TAG], state: { isCollapsed: true } }).commit();
            group = newGroupNode.ref;
        }
        return group;
    }

    async createSegmentationRepresentation3D(segmentationNode: StateObjectSelector<PluginStateObject.Volume.Data>, params: ProjectSegmentationDataParamsValues) {
        const segmentationData = segmentationNode.data as Volume;
        const segmentationId = params.segmentationId;
        const segmentation = Volume.Segmentation.get(segmentationData);
        const segmentIds: number[] = Array.from(segmentation?.segments.keys() ?? []);
        const segmentationRepresentation3D = await this.entryData.newUpdate().to(segmentationNode)
            .apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.entryData.plugin, segmentationData, {
                type: 'segment',
                typeParams: { tryUseGpu: VolsegGlobalStateData.getGlobalState(this.entryData.plugin)?.tryUseGpu },
                color: 'volume-segment',
                colorParams: { palette: this.createPalette(segmentIds) },
            }), { tags: [SEGMENT_VISUAL_TAG, segmentationId] }).commit();
    }

    // creates colors for lattice segments
    private createPalette(segmentIds: number[]) {
        const colorMap = new Map<number, Color>();
        // for (const segment of this.entryData.metadata.allSegments) {
        if (this.entryData.metadata.allAnnotations) {
            for (const annotation of this.entryData.metadata.allAnnotations) {
                if (annotation.color) {
                    const color = Color.fromNormalizedArray(annotation.color, 0);
                    colorMap.set(annotation.segment_id, color);
                }
            }
            if (colorMap.size === 0) return undefined;
            for (const segid of segmentIds) {
                colorMap.get(segid);
            }
            const colors = segmentIds.map(segid => colorMap.get(segid) ?? DEFAULT_SEGMENT_COLOR);
            return { name: 'colors' as const, params: { list: { kind: 'set' as const, colors: colors } } };
        }
    }

    async updateOpacity(opacity: number) {
        const reprs = this.entryData.findNodesByTags(SEGMENT_VISUAL_TAG);
        const update = this.entryData.newUpdate();
        for (const s of reprs) {
            update.to(s).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.alpha = opacity; });
        }
        return await update.commit();
    }
    private makeLoci(segments: number[], segmentationId: string) {
        const vis = this.entryData.findNodesByTags(SEGMENT_VISUAL_TAG, segmentationId)[0];
        if (!vis) return undefined;
        const repr = vis.obj?.data.repr;
        const wholeLoci = repr.getAllLoci()[0];
        if (!wholeLoci || !Volume.Segment.isLoci(wholeLoci)) return undefined;
        return { loci: Volume.Segment.Loci(wholeLoci.volume, segments), repr: repr };
    }
    async highlightSegment(segmentId: number, segmentationId: string) {
        const segmentLoci = this.makeLoci([segmentId], segmentationId);
        if (!segmentLoci) return;
        this.entryData.plugin.managers.interactivity.lociHighlights.highlight(segmentLoci, false);
    }
    async selectSegment(segmentId?: number, segmentationId?: string) {
        if (segmentId === undefined || segmentId < 0 || !segmentationId) return;
        const segmentLoci = this.makeLoci([segmentId], segmentationId);
        if (!segmentLoci) return;
        this.entryData.plugin.managers.interactivity.lociSelects.select(segmentLoci, false);
    }

    /** Make visible the specified set of lattice segments */
    async showSegments(segmentIds: number[], segmentationId: string) {
        const repr = this.entryData.findNodesByTags(SEGMENT_VISUAL_TAG, segmentationId)[0];
        if (!repr) return;
        const selectedSegmentKey = this.entryData.currentState.value.selectedSegment;
        const parsedSelectedSegmentKey = parseSegmentKey(selectedSegmentKey);
        const selectedSegmentId = parsedSelectedSegmentKey.segmentId;
        const selectedSegmentSegmentationId = parsedSelectedSegmentKey.segmentationId;
        const mustReselect = segmentIds.includes(selectedSegmentId) && !repr.params?.values.type.params.segments.includes(selectedSegmentId);
        const update = this.entryData.newUpdate();
        update.to(repr).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.segments = segmentIds; });
        await update.commit();
        if (mustReselect) {
            await this.selectSegment(selectedSegmentId, selectedSegmentSegmentationId);
        }
    }

    async setTryUseGpu(tryUseGpu: boolean) {
        const visuals = this.entryData.findNodesByTags(SEGMENT_VISUAL_TAG);
        for (const visual of visuals) {
            const oldParams: VolumeVisualParams = visual.transform.params;
            if (oldParams.type.params.tryUseGpu === !tryUseGpu) {
                const newParams = { ...oldParams, type: { ...oldParams.type, params: { ...oldParams.type.params, tryUseGpu: tryUseGpu } } };
                const update = this.entryData.newUpdate().to(visual.transform.ref).update(newParams);
                await PluginCommands.State.Update(this.entryData.plugin, { state: this.entryData.plugin.state.data, tree: update, options: { doNotUpdateCurrent: true } });
            }
        }
    }
}