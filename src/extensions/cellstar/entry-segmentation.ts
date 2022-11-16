import { BehaviorSubject } from 'rxjs';
import { Volume } from '../../mol-model/volume';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Download } from '../../mol-plugin-state/transforms/data';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { Color } from '../../mol-util/color';

import { Segment } from './cellstar-api/data';
import { CellStarEntryData, MAX_VOXELS } from './entry-root';
import { CreateVolume, NodeManager } from './helpers';
import { LatticeSegmentation } from './lattice-segmentation';


export class CellStarSegmentationData {
    private entryData: CellStarEntryData;

    private segmentation?: LatticeSegmentation;
    private segmentationNodeMgr = new NodeManager();
    currentSegment = new BehaviorSubject<Segment | undefined>(undefined);

    constructor(rootData: CellStarEntryData) {
        this.entryData = rootData;
    }

    async showLatticeSegmentation() {
        const hasLattices = this.entryData.metadata.grid.segmentation_lattices.segmentation_lattice_ids.length > 0;
        if (hasLattices) {
            const url = this.entryData.api.latticeUrl(this.entryData.source, this.entryData.entryId, 0, null, MAX_VOXELS);
            const data = await this.entryData.newUpdate().apply(Download, { url, isBinary: true, label: `Segmentation Data: ${url}` }).commit();
            const cif = await this.entryData.plugin.build().to(data).apply(StateTransforms.Data.ParseCif).commit();
            const latticeBlock = cif.data!.blocks.find(b => b.header === 'SEGMENTATION_DATA');
            if (latticeBlock) {
                this.segmentation = await LatticeSegmentation.fromCifBlock(latticeBlock);
                await this.showSegments(this.entryData.metadata.annotation.segment_list);
            } else {
                console.log('WARNING: Block SEGMENTATION_DATA is missing. Not showing segmentations.');
            }
        }
    }

    /** Make visible the specified set of lattice segments */
    async showSegments(segments: Segment[]) {
        this.currentSegment.next(segments.length === 1 ? segments[0] : undefined);

        const update = this.entryData.newUpdate();
        const group = await this.entryData.groupNodeMgr.showNode('Segmentation', () => update.apply(CreateGroup, {label: 'Segmentation'}).selector, false);

        this.segmentationNodeMgr.hideAllNodes();

        for (const seg of segments) {
            this.segmentationNodeMgr.showNode(seg.id.toString(), () => {
                const volume = this.segmentation?.createSegment(seg.id);
                Volume.PickingGranularity.set(volume!, 'volume');
                const volumeNode = update.to(group).apply(CreateVolume, { volume, label: `Segment ${seg.id}`, description: seg.biological_annotation?.name }, { state: { isCollapsed: true } });

                volumeNode.apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.entryData.plugin, volume, {
                    type: 'isosurface',
                    typeParams: { alpha: 1, isoValue: Volume.IsoValue.absolute(0.95) },
                    color: 'uniform',
                    colorParams: { value: Color.fromNormalizedArray(seg.colour, 0) }
                }));
                return volumeNode.selector;
            });
        }
        await update.commit();
    }
}