/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginContext } from '../../../../mol-plugin/context';
import { Color } from '../../../../mol-util/color';
import { parseCVSXJSON } from '../../common';
import { CVSXGeometricSegmentationData, CVSXMeshSegmentationData, QueryArgs } from '../../cvsx-extension/data';
import { objectToArray } from '../helpers';
import { DescriptionData, GridMetadata, Metadata, ParsedSegmentKey, ShapePrimitiveData } from './data';


export function compareTwoObjects(object1: any, object2: any): boolean {
    if (!object1 || !object2) return false;
    let compareRes = true;
    if (Object.keys(object1).length === Object.keys(object2).length) {
        Object.keys(object1).forEach(key => {
            if (object1[key] !== object2[key]) {
                compareRes = false;
            }
        });
    } else {
        compareRes = false;
    }
    return compareRes;
}

export async function getCVSXGeometricSegmentationDataFromRaw(rawData: [string, Uint8Array][], metadata: GridMetadata, query: QueryArgs, plugin: PluginContext) {
    const hasGeometricSegmentation = metadata.geometric_segmentation;
    const data: CVSXGeometricSegmentationData[] = [];
    if (hasGeometricSegmentation && hasGeometricSegmentation.segmentation_ids.length > 0) {
        let segmentationIds = hasGeometricSegmentation.segmentation_ids;
        if (query.segmentation_id) {
            segmentationIds = [query.segmentation_id];
        }
        for (const segmentationId of segmentationIds) {
            const timeInfo = hasGeometricSegmentation.time_info[segmentationId];
            let timeframes = Array.from({ length: timeInfo.end - timeInfo.start + 1 }, (_, i) => i + timeInfo.start);
            if (query.time) {
                timeframes = [query.time];
            }
            for (const timeframe of timeframes) {
                const gsData = rawData.find(r => {
                    const filename = r[0].split('.')[0];
                    const targetSegmentationId = filename.split('_')[1];
                    const timeframeIndex = parseInt(filename.split('_')[2]);
                    if (segmentationId === targetSegmentationId && timeframe === timeframeIndex) {
                        return true;
                    }
                });
                debugger;
                const parsedGsData: ShapePrimitiveData = await parseCVSXJSON(gsData!, plugin);
                const d: CVSXGeometricSegmentationData = {
                    segmentationId: segmentationId,
                    timeframeIndex: timeframe,
                    data: parsedGsData
                };
                data.push(d);
            }
        }
        return data;
    }
}
// TODO: this should consider only the files
export function getCVSXMeshSegmentationDataFromRaw(rawData: [string, Uint8Array][], metadata: GridMetadata, query: QueryArgs) {
    const hasMeshSegmentations = metadata.segmentation_meshes;
    const data: CVSXMeshSegmentationData[] = [];
    if (hasMeshSegmentations && hasMeshSegmentations.segmentation_ids.length > 0) {
        let segmentationIds = hasMeshSegmentations.segmentation_ids;
        if (query.segmentation_id) {
            segmentationIds = [query.segmentation_id];
        }
        for (const segmentationId of segmentationIds) {
            const timeInfo = hasMeshSegmentations.time_info[segmentationId];
            const timeframes = Array.from({ length: timeInfo.end - timeInfo.start + 1 }, (_, i) => i + timeInfo.start);
            for (const timeframe of timeframes) {
                const targetSegments = rawData.filter(r => {
                    const filename = r[0].split('.')[0];
                    // const segmentId = filename.split('_')[1];
                    const targetSegmentationId = filename.split('_')[2];
                    const timeframeIndex = parseInt(filename.split('_')[3]);
                    if (segmentationId === targetSegmentationId && timeframe === timeframeIndex) {
                        return true;
                    }
                });
                console.log(`Target Segments for ${timeframe} and ${segmentationId}`);
                console.log(targetSegments);
                const d: CVSXMeshSegmentationData = {
                    segmentationId: segmentationId,
                    timeframeIndex: timeframe,
                    data: targetSegments
                };
                data.push(d);
            }
        }
        return data;
    } else {
        return false;
    }
}

export class MetadataWrapper {
    raw: Metadata;
    private segmentMap?: Map<string, DescriptionData[]>;

    constructor(rawMetadata: Metadata) {
        this.raw = rawMetadata;
    }

    hasLatticeSegmentations() {
        const grid = this.raw.grid;
        if (grid.segmentation_lattices && grid.segmentation_lattices.segmentation_ids.length > 0) {
            return grid.segmentation_lattices;
        }
        return false;
    };

    hasMeshSegmentations() {
        const grid = this.raw.grid;
        if (grid.segmentation_meshes && grid.segmentation_meshes.segmentation_ids.length > 0) {
            return grid.segmentation_meshes;
        }
        return false;
    };

    hasGeometricSegmentations() {
        const grid = this.raw.grid;
        if (grid.geometric_segmentation && grid.geometric_segmentation.segmentation_ids.length > 0) {
            return grid.geometric_segmentation;
        }
        return false;
    };

    hasSegmentations() {
        const grid = this.raw.grid;
        if (grid.geometric_segmentation && grid.geometric_segmentation.segmentation_ids.length > 0) {
            return true;
        }
        if (grid.segmentation_lattices && grid.segmentation_lattices.segmentation_ids.length > 0) {
            return true;
        }
        if (grid.segmentation_meshes && grid.segmentation_meshes.segmentation_ids.length > 0) {
            return true;
        }
        return false;
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
            return (objectToArray(descriptions) as DescriptionData[]);
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

    getAllDescriptionsBasedOnMetadata(metadata: { [key: string]: any }) {
        const a = this.allDescriptions;
        const filtered = a.filter(i => compareTwoObjects(i.metadata, metadata));
        return filtered;
    }

    filterDescriptionsBasedOnKeywordInMetadata(d: DescriptionData[], keyword: string) {
        // for each description
        // loop over metadata fields
        // compare value of each field with keyword
        if (keyword === '') return d;
        const filtered = d.filter(i => i.metadata && Object.values(i.metadata).includes(keyword));
        // suppose no segments
        return filtered;
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

export function instanceOfShapePrimitiveData(object: any): object is ShapePrimitiveData {
    debugger;
    return 'shape_primitive_list' in object;
}

export function createSegmentKey(segmentId: number, segmentationId: string, kind: 'lattice' | 'mesh' | 'primitive') {
    return `${kind}:${segmentationId}:${segmentId}`;
}

export function parseSegmentKey(segmentKey: string) {
    const kind = segmentKey.split(':')[0] as 'lattice' | 'mesh' | 'primitive';
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
            { id: description.target_id!.segment_id, label: description?.name && !description?.is_hidden ? `<b>${description?.name}</b>` : 'No label provided or description is hidden' }
        )
    );
}
