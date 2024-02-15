/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Color } from '../../../mol-util/color';
import { DescriptionData, Metadata, ParsedSegmentKey } from './data';

export class MetadataWrapper {
    raw: Metadata;
    private segmentMap?: Map<string, DescriptionData[]>;

    constructor(rawMetadata: Metadata) {
        this.raw = rawMetadata;
    }

    removeDescription(id: string) {
        const descriptions = this.raw.annotation?.descriptions;
        if (this.raw.annotation && descriptions) {
            // const removed = descriptions.filter(d => d.i)
            delete descriptions[id];
            this.raw.annotation.descriptions = descriptions;
            console.log(this.raw.annotation.descriptions);
        }
    }

    removeSegmentAnnotation(id: string) {
        const segmentAnnotations = this.allAnnotations;
        const filtered = segmentAnnotations.filter(a => a.id !== id);
        this.raw.annotation!.annotations = filtered;
    }

    get allDescriptions() {
        const descriptions = this.raw.annotation?.descriptions;
        if (descriptions) {
            const d = [];
            const arr = Object.entries(descriptions);
            for (const obj of arr) {
                d.push(obj[1]);
            };
            return d;
        } else {
            return [];
        }
    }

    getAllDescriptionsForSegmentationAndTimeframe(segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive', timeframeIndex: number) {
        const allDescriptions = this.allDescriptions;
        const allDescriptonsForSegmentationId = allDescriptions.filter(d =>
            d.target_id && d.target_id.segmentation_id === segmentationId && d.target_kind === kind
        );
        const allDescriptonsForSegmentationIdAndTimeframe = allDescriptonsForSegmentationId.filter(a => {
            if (a.hasOwnProperty('time') && Number.isFinite(a.time)) {
                if (timeframeIndex === a.time) {
                    return a;
                }
            } else if (a.hasOwnProperty('time') && Array.isArray(a.time) && a.time.every(i => Number.isFinite(i))) {
                if ((a.time as number[]).includes(timeframeIndex)) {
                    return a;
                }
            } else {
                return a;
            }
        });

        return allDescriptonsForSegmentationIdAndTimeframe;
    }

    get allAnnotations() {
        const annotations = this.raw.annotation?.annotations;
        if (annotations) return annotations; else return [];
    }
    getAllAnnotationsForTimeframe(timeframeIndex: number) {
        const allSegmentsForTimeframe = this.allAnnotations.filter(a => {
            if (a.hasOwnProperty('time') && Number.isFinite(a.time)) {
                if (timeframeIndex === a.time) {
                    return a;
                }
            } else if (a.hasOwnProperty('time') && Array.isArray(a.time) && a.time.every(i => Number.isFinite(i))) {
                if ((a.time as number[]).includes(timeframeIndex)) {
                    return a;
                }
            } else {
                return a;
            }
        });
        return allSegmentsForTimeframe;
    }

    get allSegmentIds() {
        // return this.allSegments.map(segment => segment.id);
        // TODO: Segment ids should be unique
        return this.allAnnotations?.map(annotation => annotation.segment_id);
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

    getSegmentAnnotation(segmentId: number, segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive') {
        const allAnnotations = this.allAnnotations;
        return allAnnotations.find(a => a.segment_id === segmentId && a.segment_kind === kind && a.segmentation_id === segmentationId);
    }

    getSegmentDescription(segmentId: number, segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive'): DescriptionData[] | undefined {
        // NOTE: for now assumes single lattice, single mesh set etc.
        const segmentKey = createSegmentKey(segmentId, segmentationId, kind);
        if (this.allDescriptions) {
            // create segment map
            if (!this.segmentMap) {
                this.segmentMap = new Map<string, DescriptionData[]>;
                for (const description of this.allDescriptions) {
                    if (description.target_id && description.target_kind !== 'entry') {
                        const _kind = description.target_kind;
                        const _segmentationId = description.target_id.segmentation_id;
                        const _segmentId = description.target_id.segment_id;
                        const _segmentKey = createSegmentKey(_segmentId, _segmentationId, _kind);
                        const existingDescriptions = this.segmentMap.get(_segmentKey);
                        if (existingDescriptions) {
                            this.segmentMap.set(_segmentKey, [...existingDescriptions, description]);
                        } else {
                            this.segmentMap.set(_segmentKey, [description]);
                        }
                    }
                }
            }
            return this.segmentMap.get(segmentKey);
        }
    }

    // TODO: apparently for meshes
    // getSegmentColor(segmentId: number): Color | undefined {
    //     const colorArray = this.getSegmentDescription(segmentId)?.color;
    //     return colorArray ? Color.fromNormalizedArray(colorArray, 0) : undefined;
    // }

    /** Get the list of detail levels available for the given mesh segment. */
    getMeshDetailLevels(segmentationId: string, timeframe: number, segmentId: number): number[] {
        if (!this.raw.grid.segmentation_meshes) return [];
        const meshComponentNumbers = this.raw.grid.segmentation_meshes.segmentation_metadata[segmentationId].mesh_timeframes[timeframe];
        // const segmentIds = segmentationSetMetadata.

        const segmentIds = meshComponentNumbers.segment_ids;
        if (!segmentIds) return [];
        const details = segmentIds[segmentId].detail_lvls;
        return Object.keys(details).map(s => parseInt(s));
        // const segmentIds = this.raw.grid.segmentation_meshes.mesh_component_numbers.segment_ids;
        // if (!segmentIds) return [];
        // const details = segmentIds[segmentId].detail_lvls;
        // return Object.keys(details).map(s => parseInt(s));
    }

    /** Get the worst available detail level that is not worse than preferredDetail.
     * If preferredDetail is null, get the worst detail level overall.
     * (worse = greater number) */
    getSufficientMeshDetail(segmentationId: string, timeframe: number, segmentId: number, preferredDetail: number | null) {
        let availDetails = this.getMeshDetailLevels(segmentationId, timeframe, segmentId);
        if (preferredDetail !== null) {
            availDetails = availDetails.filter(det => det <= preferredDetail);
        }
        return Math.max(...availDetails);
    }

    /** IDs of all segments available as meshes */
    // make it work for multiple mesh sets?
    getMeshSegmentIdsForSegmentationIdAndTimeframe(segmentationId: string, timeframe: number) {
        // if (!this.raw.grid.segmentation_meshes) return;
        const meshComponentNumbers = this.raw.grid.segmentation_meshes!.segmentation_metadata[segmentationId].mesh_timeframes[timeframe];
        // const segmentIds = segmentationSetMetadata.

        const segmentIds = meshComponentNumbers.segment_ids;
        if (!segmentIds) return [];
        return Object.keys(segmentIds).map(s => parseInt(s));
    }

    getAllSegmentAnotationsForSegmentationAndTimeframe(segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive', timeframeIndex: number) {
        const allAnnotations = this.allAnnotations;
        const allAnnotationsForSegmentationId = allAnnotations.filter(d =>
            d.segmentation_id === segmentationId && d.segment_kind === kind
        );
        const allAnnotationsForSegmentationIdAndTimeframe = allAnnotationsForSegmentationId.filter(a => {
            if (a.hasOwnProperty('time') && Number.isFinite(a.time)) {
                if (timeframeIndex === a.time) {
                    return a;
                }
            } else if (a.hasOwnProperty('time') && Array.isArray(a.time) && a.time.every(i => Number.isFinite(i))) {
                if ((a.time as number[]).includes(timeframeIndex)) {
                    return a;
                }
            } else {
                return a;
            }
        });

        return allAnnotationsForSegmentationIdAndTimeframe;
    }

    // get meshSegmentIds() {
    //     const segmentIds = this.raw.grid.segmentation_meshes.mesh_component_numbers.segment_ids;
    //     if (!segmentIds) return [];
    //     return Object.keys(segmentIds).map(s => parseInt(s));
    // }

    get gridTotalVolume() {
        const sortedResolutions = this.raw.grid.volumes.volume_sampling_info.spatial_downsampling_levels.sort((n1, n2) => n1 - n2);
        const firstAvailableResolution = sortedResolutions[0];
        console.log('firstAvailableResolution', firstAvailableResolution);
        const [vx, vy, vz] = this.raw.grid.volumes.volume_sampling_info.boxes[firstAvailableResolution].voxel_size;
        const [gx, gy, gz] = this.raw.grid.volumes.volume_sampling_info.boxes[firstAvailableResolution].grid_dimensions;
        return vx * vy * vz * gx * gy * gz;
    }

}

export function createSegmentKey(segmentId: number, segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive') {
    return `${kind}:${segmentationId}:${segmentId}`;
}

export function parseSegmentKey(segmentKey: string) {
    const kind: 'lattice' | 'mesh' | 'primitive' = segmentKey.split(':')[0];
    const segmentationId = segmentKey.split(':')[1];
    const segmentId = segmentKey.split(':')[2];
    const parsedSegmentKey: ParsedSegmentKey = {
        kind: kind,
        segmentationId: segmentationId,
        segmentId: Number(segmentId)
    };
    return parsedSegmentKey;
}
export function getSegmentLabelsFromDescriptions(descriptions: DescriptionData[]) {
    return descriptions.map(
        description => (
            { id: description.target_id!.segment_id, label: description?.name ? description?.name : 'No label provided' }
        )
    );
}
