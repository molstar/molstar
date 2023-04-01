/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Color } from '../../../mol-util/color';
import { Metadata, Segment } from './data';


export class MetadataWrapper {
    raw: Metadata;
    private segmentMap?: { [id: number]: Segment };

    constructor(rawMetadata: Metadata) {
        this.raw = rawMetadata;
    }

    get allSegments() {
        return this.raw.annotation?.segment_list ?? [];
    }

    get allSegmentIds() {
        return this.allSegments.map(segment => segment.id);
    }

    getSegment(segmentId: number): Segment | undefined {
        if (!this.segmentMap) {
            this.segmentMap = {};
            for (const segment of this.allSegments) {
                this.segmentMap[segment.id] = segment;
            }
        }
        return this.segmentMap[segmentId];
    }

    getSegmentColor(segmentId: number): Color | undefined {
        const colorArray = this.getSegment(segmentId)?.colour;
        return colorArray ? Color.fromNormalizedArray(colorArray, 0) : undefined;
    }

    /** Get the list of detail levels available for the given mesh segment. */
    getMeshDetailLevels(segmentId: number): number[] {
        const segmentIds = this.raw.grid.segmentation_meshes.mesh_component_numbers.segment_ids;
        if (!segmentIds) return [];
        const details = segmentIds[segmentId].detail_lvls;
        return Object.keys(details).map(s => parseInt(s));
    }

    /** Get the worst available detail level that is not worse than preferredDetail.
     * If preferredDetail is null, get the worst detail level overall.
     * (worse = greater number) */
    getSufficientMeshDetail(segmentId: number, preferredDetail: number | null) {
        let availDetails = this.getMeshDetailLevels(segmentId);
        if (preferredDetail !== null) {
            availDetails = availDetails.filter(det => det <= preferredDetail);
        }
        return Math.max(...availDetails);
    }

    /** IDs of all segments available as meshes */
    get meshSegmentIds() {
        const segmentIds = this.raw.grid.segmentation_meshes.mesh_component_numbers.segment_ids;
        if (!segmentIds) return [];
        return Object.keys(segmentIds).map(s => parseInt(s));
    }

    get gridTotalVolume() {
        const [vx, vy, vz] = this.raw.grid.volumes.voxel_size[1];
        const [gx, gy, gz] = this.raw.grid.volumes.grid_dimensions;
        return vx * vy * vz * gx * gy * gz;
    }

}
