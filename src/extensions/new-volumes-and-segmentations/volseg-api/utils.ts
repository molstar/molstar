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
        let allSegments: Segment[] = [];
        const segmentationLattices = this.raw.annotation?.segmentation_lattices;
        if (segmentationLattices) {
            for (const segmentationLatticeInfo of segmentationLattices) {
                allSegments = allSegments.concat(segmentationLatticeInfo.segment_list);
            };
        }
        return allSegments;
    }

    get allSegmentIds() {
        return this.allSegments.map(segment => segment.id);
    }

    get channelAnnotations() {
        return this.raw.annotation?.volume_channels_annotations;
    }

    getVolumeChannelColor(channel_id: string) {
        if (!this.channelAnnotations) {
            return Color(0x121212);
        }
        const channelColorArray = this.channelAnnotations.filter(i => i.channel_id === channel_id)[0]?.color;
        if (channelColorArray) {
            const color = Color.fromNormalizedArray(channelColorArray, 0);
            return color;
        } else {
            return Color(0x121212);
        }
    }

    getVolumeChannelLabel(channel_id: string) {
        if (!this.channelAnnotations) {
            return null;
        }
        const volumeChannelLabel = this.channelAnnotations
            .filter(i => i.channel_id === channel_id)[0]?.label;

        if (volumeChannelLabel) {
            return volumeChannelLabel;
        } else {
            return null;
        };

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
        const colorArray = this.getSegment(segmentId)?.color;
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
        const [vx, vy, vz] = this.raw.grid.volumes.volume_sampling_info.boxes[1].voxel_size;
        const [gx, gy, gz] = this.raw.grid.volumes.volume_sampling_info.boxes[1].grid_dimensions;
        return vx * vy * vz * gx * gy * gz;
    }

}
