import { CIF } from '../../mol-io/reader/cif';
import { Volume } from '../../mol-model/volume';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { Asset } from '../../mol-util/assets';
import { Color } from '../../mol-util/color';
import { PluginCommands, PluginStateObject } from '../meshes/molstar-lib-imports';

import { Segment } from './cellstar-api/data';
import { CellStarEntryData, MAX_VOXELS } from './entry-root';
import { CreateVolume, NodeManager } from './helpers';
import { LatticeSegmentation } from './lattice-segmentation';


const LATTICE_SEGMENT_TAG = 'lattice-segment';


export class CellStarLatticeSegmentationData {
    private entryData: CellStarEntryData;

    private segmentation?: LatticeSegmentation;
    private segmentationNodeMgr = new NodeManager();


    constructor(rootData: CellStarEntryData) {
        this.entryData = rootData;
    }

    async showSegmentation() {
        const hasLattices = this.entryData.metadata.grid.segmentation_lattices.segmentation_lattice_ids.length > 0;
        if (hasLattices) {
            const url = this.entryData.api.latticeUrl(this.entryData.source, this.entryData.entryId, 0, null, MAX_VOXELS);

            const urlAsset = Asset.getUrlAsset(this.entryData.plugin.managers.asset, url);
            const asset = await this.entryData.plugin.runTask(this.entryData.plugin.managers.asset.resolve(urlAsset, 'binary'));
            const parsed = await this.entryData.plugin.runTask(CIF.parseBinary(asset.data));
            if (parsed.isError) {
                throw new Error(`Failed parsing CIF file from ${url}`);
            }
            const latticeBlock = parsed.result.blocks.find(b => b.header === 'SEGMENTATION_DATA');
            if (latticeBlock) {
                this.segmentation = await LatticeSegmentation.fromCifBlock(latticeBlock);
                await this.showSegments(this.entryData.metadata.annotation?.segment_list ?? []);
            } else {
                console.log('WARNING: Block SEGMENTATION_DATA is missing. Not showing segmentations.');
            }
        }
    }

    updateOpacity(opacity: number) {
        const group = this.entryData.groupNodeMgr.getNode('LatticeSegmentation');
        if (!group) return;

        const segs = this.entryData.plugin.state.data.selectQ(q => q.byRef(group.ref).subtree().withTag(LATTICE_SEGMENT_TAG));
        const update = this.entryData.newUpdate();
        for (const s of segs) {
            update.to(s).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.alpha = opacity; })
        }
        return update.commit();
    }

    highlightSegment(segment?: Segment) {
        PluginCommands.Interactivity.ClearHighlights(this.entryData.plugin);
        if (!segment) return;
        const node = this.segmentationNodeMgr.getNode(segment.id.toString());
        if (!node) return;
        const vis = this.entryData.plugin.state.data.selectQ(q => q.byRef(node.ref).subtree().ofType(PluginStateObject.Volume.Representation3D))[0];
        if (!vis) return;
        PluginCommands.Interactivity.Object.Highlight(this.entryData.plugin, { state: this.entryData.plugin.state.data, ref: vis.transform.ref });
    }


    /** Make visible the specified set of lattice segments */
    async showSegments(segments: Segment[], options?: { opacity?: number }) {
        this.segmentationNodeMgr.hideAllNodes();

        segments = segments.filter(seg => this.segmentation?.hasSegment(seg.id));
        if (segments.length == 0) return;

        const group = await this.entryData.groupNodeMgr.showNode('LatticeSegmentation', async () => await this.entryData.newUpdate().apply(CreateGroup, { label: 'Segmentation', description: 'Lattice' }, { state: { isCollapsed: true } }).commit(), false)

        const update = this.entryData.newUpdate();
        for (const seg of segments) {
            this.segmentationNodeMgr.showNode(seg.id.toString(), () => {
                const volume = this.segmentation?.createSegment(seg, { ownerId: this.entryData.entryRoot?.ref, segment: seg });
                Volume.PickingGranularity.set(volume!, 'volume');
                const volumeNode = update.to(group).apply(CreateVolume, { volume, label: `Segment ${seg.id}`, description: seg.biological_annotation?.name }, { state: { isCollapsed: true } });

                volumeNode.apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.entryData.plugin, volume, {
                    type: 'isosurface',
                    typeParams: { alpha: options?.opacity ?? 1, isoValue: Volume.IsoValue.absolute(0.95) },
                    color: 'uniform',
                    colorParams: { value: Color.fromNormalizedArray(seg.colour, 0) }
                }), { tags: [LATTICE_SEGMENT_TAG] });
                return volumeNode.selector;
            }, undefined);
        }
        await update.commit();
    }
}