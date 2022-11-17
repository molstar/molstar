import { BehaviorSubject } from 'rxjs';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { Color } from '../../mol-util/color';
import { ColorNames } from '../../mol-util/color/names';
import { createMeshFromUrl } from '../meshes/examples';
import { BACKGROUND_SEGMENT_VOLUME_THRESHOLD } from '../meshes/mesh-streaming/behavior';

import { Segment } from './cellstar-api/data';
import { CellStarEntryData } from './entry-root';
import { MetadataUtils, NodeManager } from './helpers';


const DEFAULT_MESH_DETAIL: number | null = 5;  // null means worst


export class CellStarMeshSegmentationData {
    private entryData: CellStarEntryData;

    private segmentationNodeMgr = new NodeManager();
    currentSegment = new BehaviorSubject<Segment | undefined>(undefined);

    constructor(rootData: CellStarEntryData) {
        this.entryData = rootData;
    }

    async showSegmentation() {
        const hasMeshes = this.entryData.metadata.grid.segmentation_meshes.mesh_component_numbers.segment_ids !== undefined;
        if (hasMeshes) {
            await this.showSegments(this.entryData.metadata.annotation.segment_list);
        }
    }

    /** Make visible the specified set of mesh segments */
    async showSegments(segments: Segment[]) {
        this.currentSegment.next(segments.length === 1 ? segments[0] : undefined);
        this.segmentationNodeMgr.hideAllNodes();

        const meshSegments = this.entryData.metadata.grid.segmentation_meshes.mesh_component_numbers.segment_ids ?? {};

        segments = segments.filter(seg => meshSegments[seg.id] !== undefined);
        if (segments.length == 0) return;

        const group = await this.entryData.groupNodeMgr.showNode('MeshSegmentation', async () => await this.entryData.newUpdate().apply(CreateGroup, { label: 'Segmentation', description: 'Mesh' }, { state: { isCollapsed: true } }).commit(), false)

        const [vx, vy, vz] = this.entryData.metadata.grid.volumes.voxel_size[1];
        const [gx, gy, gz] = this.entryData.metadata.grid.volumes.grid_dimensions;
        const totalVolume = vx * vy * vz * gx * gy * gz;

        for (const seg of segments) {
            this.segmentationNodeMgr.showNode(seg.id.toString(), async () => {
                const detail = MetadataUtils.getSufficientDetail(this.entryData.metadata!, seg.id, DEFAULT_MESH_DETAIL);
                const color = seg.colour.length >= 3 ? Color.fromNormalizedArray(seg.colour, 0) : ColorNames.gray;
                return await createMeshFromUrl(this.entryData.plugin, this.entryData.api.meshUrl_Bcif(this.entryData.source, this.entryData.entryId, seg.id, detail), seg.id, detail, 
                    true, false, color, group, BACKGROUND_SEGMENT_VOLUME_THRESHOLD * totalVolume);
            });
        }
    }
}