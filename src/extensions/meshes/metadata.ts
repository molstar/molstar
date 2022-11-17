import * as MS from './molstar-lib-imports';


// TODO unify with Metadata in CellStar

export interface Metadata {
    grid: {
        general: {
            details: string,
        },
        volumes: Volumes,
        segmentation_lattices: SegmentationLattices,
        segmentation_meshes: SegmentationMeshes,
    },
    annotation: Annotation,
}

export interface Volumes {
    volume_downsamplings: number[],
    voxel_size: { [downsampling: number]: Vector3 },
    origin: Vector3,
    grid_dimensions: Vector3,
    sampled_grid_dimensions: { [downsampling: number]: Vector3 },
    mean: { [downsampling: number]: number },
    std: { [downsampling: number]: number },
    min: { [downsampling: number]: number },
    max: { [downsampling: number]: number },
    volume_force_dtype: string,
}

export interface SegmentationLattices {
    segmentation_lattice_ids: number[],
    segmentation_downsamplings: { [lattice: number]: number[] },
}

export interface Annotation {
    name: string,
    details: string,
    segment_list: Segment[],
}

export interface Segment {
    id: number,
    colour: number[],
    biological_annotation: BiologicalAnnotation,
}

export interface BiologicalAnnotation {
    name: string,
    external_references: { id: number, resource: string, accession: string, label: string, description: string }[]
}

export interface SegmentationMeshes {
    mesh_component_numbers: {
        segment_ids?: {
            [segId: number]: {
                detail_lvls: {
                    [detail: number]: {
                        mesh_ids: {
                            [meshId: number]: {
                                num_triangles: number,
                                num_vertices: number
                            }
                        }
                    }
                }
            }
        }
    }
    detail_lvl_to_fraction: {
        [lvl: number]: number
    }
}

type Vector3 = [number, number, number];



export namespace Metadata {
    export function meshSegments(metadata: Metadata): number[] {
        const segmentIds = metadata.grid.segmentation_meshes.mesh_component_numbers.segment_ids;
        if (segmentIds === undefined) return [];
        return Object.keys(segmentIds).map(s => parseInt(s));
    }
    export function meshSegmentDetails(metadata: Metadata, segmentId: number): number[] {
        const segmentIds = metadata.grid.segmentation_meshes.mesh_component_numbers.segment_ids;
        if (segmentIds === undefined) return [];
        const details = segmentIds[segmentId].detail_lvls;
        return Object.keys(details).map(s => parseInt(s));
    }
    /** Get the worst available detail level that is not worse than preferredDetail.
     * If preferredDetail is null, get the worst detail level overall.
     * (worse = greater number) */
    export function getSufficientDetail(metadata: Metadata, segmentId: number, preferredDetail: number | null) {
        let availDetails = meshSegmentDetails(metadata, segmentId);
        if (preferredDetail !== null) {
            availDetails = availDetails.filter(det => det <= preferredDetail);
        }
        return Math.max(...availDetails);
    }
    export function annotationsBySegment(metadata: Metadata): { [id: number]: Segment } {
        const result: { [id: number]: Segment } = {};
        for (const segment of metadata.annotation.segment_list) {
            if (segment.id in result) {
                throw new Error(`Duplicate segment annotation for segment ${segment.id}`);
            }
            result[segment.id] = segment;
        }
        return result;
    }
    export function dropSegments(metadata: Metadata, segments: number[]): void {
        if (metadata.grid.segmentation_meshes.mesh_component_numbers.segment_ids === undefined) return;
        const dropSet = new Set(segments);
        metadata.annotation.segment_list = metadata.annotation.segment_list.filter(seg => !dropSet.has(seg.id));
        for (const seg of segments) {
            delete metadata.grid.segmentation_meshes.mesh_component_numbers.segment_ids[seg];
        }
    }
    export function namesAndColorsBySegment(metadata: Metadata) {
        const result: { [id: number]: { name: string, color: MS.Color } } = {};
        for (const segment of metadata.annotation.segment_list) {
            if (segment.id in result) throw new Error(`Duplicate segment annotation for segment ${segment.id}`);
            result[segment.id] = { name: segment.biological_annotation.name, color: MS.Color.fromNormalizedArray(segment.colour, 0) };
        }
        return result;

    }
}