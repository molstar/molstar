/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

export interface Metadata {
    grid: {
        general: {
            details: string,
            source_db_name: string,
            source_db_id: string,
        },
        volumes: Volumes,
        segmentation_lattices: SegmentationLattices,
        segmentation_meshes: SegmentationMeshes,
    },
    annotation: Annotation | null,
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
    external_references: ExternalReference[]
}

export interface ExternalReference {
    id: number, resource: string, accession: string, label: string,
    description: string
}

type Vector3 = [number, number, number];
