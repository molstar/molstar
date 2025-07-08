/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

import { Segment } from './volseg-api/data';
import { BOX, VolsegEntryData, MAX_VOXELS } from './entry-root';
import { VolumeVisualParams } from './entry-volume';
import { VolsegGlobalStateData } from './global-state';
import { Interval } from '../../mol-data/int/interval';
import { SortedArray } from '../../mol-data/int';


const GROUP_TAG = 'lattice-segmentation-group';
const SEGMENT_VISUAL_TAG = 'lattice-segment-visual';

const DEFAULT_SEGMENT_COLOR = Color.fromNormalizedRgb(0.8, 0.8, 0.8);


export class VolsegLatticeSegmentationData {
    private entryData: VolsegEntryData;

    constructor(rootData: VolsegEntryData) {
        this.entryData = rootData;
    }

    async loadSegmentation() {
        const hasLattices = this.entryData.metadata.raw.grid.segmentation_lattices.segmentation_lattice_ids.length > 0;
        if (hasLattices) {
            const url = this.entryData.api.latticeUrl(this.entryData.source, this.entryData.entryId, 0, BOX, MAX_VOXELS);
            let group = this.entryData.findNodesByTags(GROUP_TAG)[0]?.transform.ref;
            if (!group) {
                const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup,
                    { label: 'Segmentation', description: 'Lattice' }, { tags: [GROUP_TAG], state: { isCollapsed: true } }).commit();
                group = newGroupNode.ref;
            }
            const segmentLabels = this.entryData.metadata.allSegments.map(seg => ({ id: seg.id, label: seg.biological_annotation.name ? `<b>${seg.biological_annotation.name}</b>` : '' }));
            const volumeNode = await this.entryData.newUpdate().to(group)
                .apply(Download, { url, isBinary: true, label: `Segmentation Data: ${url}` })
                .apply(ParseCif)
                .apply(VolumeFromSegmentationCif, { blockHeader: 'SEGMENTATION_DATA', segmentLabels: segmentLabels, ownerId: this.entryData.ref })
                .commit();
            const volumeData = volumeNode.data as Volume;
            const segmentation = Volume.Segmentation.get(volumeData);
            const segmentIds: number[] = Array.from(segmentation?.segments.keys() ?? []);
            await this.entryData.newUpdate().to(volumeNode)
                .apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.entryData.plugin, volumeData, {
                    type: 'segment',
                    typeParams: { tryUseGpu: VolsegGlobalStateData.getGlobalState(this.entryData.plugin)?.tryUseGpu },
                    color: 'volume-segment',
                    colorParams: { palette: this.createPalette(segmentIds) },
                }), { tags: [SEGMENT_VISUAL_TAG] }).commit();
        }
    }

    private createPalette(segmentIds: number[]) {
        const colorMap = new Map<number, Color>();
        for (const segment of this.entryData.metadata.allSegments) {
            const color = Color.fromNormalizedArray(segment.colour, 0);
            colorMap.set(segment.id, color);
        }
        if (colorMap.size === 0) return undefined;
        for (const segid of segmentIds) {
            colorMap.get(segid);
        }
        const colors = segmentIds.map(segid => colorMap.get(segid) ?? DEFAULT_SEGMENT_COLOR);
        return { name: 'colors' as const, params: { list: { kind: 'set' as const, colors: colors } } };
    }

    async updateOpacity(opacity: number) {
        const reprs = this.entryData.findNodesByTags(SEGMENT_VISUAL_TAG);
        const update = this.entryData.newUpdate();
        for (const s of reprs) {
            update.to(s).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.alpha = opacity; });
        }
        return await update.commit();
    }
    private makeLoci(segments: number[]) {
        const vis = this.entryData.findNodesByTags(SEGMENT_VISUAL_TAG)[0];
        if (!vis) return undefined;
        const repr = vis.obj?.data.repr;
        const wholeLoci = repr.getAllLoci()[0];
        if (!wholeLoci || !Volume.Segment.isLoci(wholeLoci)) return undefined;
        const elements = [{
            segments: SortedArray.ofUnsortedArray<Volume.SegmentIndex>(segments),
            instances: Interval.ofLength(wholeLoci.volume.instances.length as Volume.InstanceIndex)
        }];
        return { loci: Volume.Segment.Loci(wholeLoci.volume, elements), repr: repr };
    }
    async highlightSegment(segment: Segment) {
        const segmentLoci = this.makeLoci([segment.id]);
        if (!segmentLoci) return;
        this.entryData.plugin.managers.interactivity.lociHighlights.highlight(segmentLoci, false);
    }
    async selectSegment(segment?: number) {
        if (segment === undefined || segment < 0) return;
        const segmentLoci = this.makeLoci([segment]);
        if (!segmentLoci) return;
        this.entryData.plugin.managers.interactivity.lociSelects.select(segmentLoci, false);
    }

    /** Make visible the specified set of lattice segments */
    async showSegments(segments: number[]) {
        const repr = this.entryData.findNodesByTags(SEGMENT_VISUAL_TAG)[0];
        if (!repr) return;
        const selectedSegment = this.entryData.currentState.value.selectedSegment;
        const mustReselect = segments.includes(selectedSegment) && !repr.params?.values.type.params.segments.includes(selectedSegment);
        const update = this.entryData.newUpdate();
        update.to(repr).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.segments = segments; });
        await update.commit();
        if (mustReselect) {
            await this.selectSegment(this.entryData.currentState.value.selectedSegment);
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