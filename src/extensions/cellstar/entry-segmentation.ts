import { Volume } from '../../mol-model/volume';
import { SegcifProvider } from '../../mol-plugin-state/formats/volume';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Download } from '../../mol-plugin-state/transforms/data';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { StateObjectSelector } from '../../mol-state';
import { Color } from '../../mol-util/color';

import { Segment } from './cellstar-api/data';
import { BOX, CellStarEntryData, MAX_VOXELS } from './entry-root';


const GROUP_NAME = 'LatticeSegmentation';
const SEGMENT_REPR_TAG = 'lattice-segment-repr';
const DEFAULT_SEGMENT_COLOR = Color.fromNormalizedRgb(0.8, 0.8, 0.8);


export class CellStarLatticeSegmentationData {
    private entryData: CellStarEntryData;

    constructor(rootData: CellStarEntryData) {
        this.entryData = rootData;
    }

    async showSegmentation() {
        const hasLattices = this.entryData.metadata.grid.segmentation_lattices.segmentation_lattice_ids.length > 0;
        if (hasLattices) {
            const url = this.entryData.api.latticeUrl(this.entryData.source, this.entryData.entryId, 0, BOX, MAX_VOXELS);
            const group = await this.entryData.groupNodeMgr.showNode(GROUP_NAME,
                async () => await this.entryData.newUpdate().apply(CreateGroup, { label: 'Segmentation', description: 'Lattice' }, { state: { isCollapsed: false } }).commit(),
                false);
            const data = await this.entryData.newUpdate().to(group).apply(Download, { url, isBinary: true, label: `Segmentation Data: ${url}` }).commit();
            console.log(this.entryData.plugin.dataFormats.list);
            const parsed = await SegcifProvider.parse(this.entryData.plugin, data);
            const volume: StateObjectSelector<PluginStateObject.Volume.Data> = parsed.volumes?.[0];
            const volumeData = volume.cell!.obj!.data;
            console.log('showSegmentation', this.entryData.ref);
            volumeData._propertyData.ownerId = this.entryData.ref;
            const segmentation = Volume.Segmentation.get(volumeData);
            if (!segmentation) return;
            segmentation.labels = {};
            for (const segment of this.entryData.metadata.annotation?.segment_list ?? []) {
                segmentation.labels[segment.id] = segment.biological_annotation.name;
            }
            const segmentIds: number[] = Array.from(segmentation.segments.keys());
            await this.entryData.plugin.build()
                .to(volume)
                .apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.entryData.plugin, volumeData, {
                    type: 'segment',
                    color: 'volume-segment',
                    colorParams: { palette: this.createPalette(segmentIds) },
                }), { tags: [SEGMENT_REPR_TAG] }).commit();
        }
    }

    private createPalette(segmentIds: number[]) {
        const colorMap = new Map<number, Color>();
        for (const segment of this.entryData.metadata.annotation?.segment_list ?? []) {
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
        const reprs = this.entryData.plugin.state.data.selectQ(q => q.byRef(this.entryData.ref).subtree().withTag(SEGMENT_REPR_TAG));
        const update = this.entryData.newUpdate();
        for (const s of reprs) {
            update.to(s).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.alpha = opacity; });
        }
        return await update.commit();
    }
    async highlightSegment(segment: Segment) {
        const vis = this.entryData.plugin.state.data.selectQ(q => q.byRef(this.entryData.ref).subtree().withTag(SEGMENT_REPR_TAG))[0];
        const repr = vis.obj?.data.repr;
        const wholeLoci = repr.getAllLoci()[0];
        if (!wholeLoci || !Volume.Segment.isLoci(wholeLoci)) return;
        const segmentLoci = Volume.Segment.Loci(wholeLoci.volume, [segment.id]);
        this.entryData.plugin.managers.interactivity.lociHighlights.highlight({ loci: segmentLoci, repr }, false);
    }

    /** Make visible the specified set of lattice segments */
    async showSegments(segments: Segment[], options?: { opacity?: number }) {
        const reprs = this.entryData.plugin.state.data.selectQ(q => q.byRef(this.entryData.ref).subtree().withTag(SEGMENT_REPR_TAG));
        const update = this.entryData.newUpdate();
        for (const repr of reprs) {
            update.to(repr).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.segments = segments.map(seg => seg.id); });
        }
        return await update.commit();
    }
}