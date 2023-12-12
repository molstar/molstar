/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

export interface Metadata {
    grid: GridMetadata,
    annotation: AnnotationMetadata | null,
}

export interface GridMetadata {
    volumes: VolumesMetadata,
    segmentation_lattices: SegmentationLatticesMetadata,
    segmentation_meshes: MeshesMetadata
    geometric_segmentation?: GeometricSegmentation | undefined
}

export interface GeometricSegmentation {
    exists?: boolean
}

export interface VolumesMetadata {
    channel_ids: number[]
    time_info: TimeInfo
    volume_sampling_info: VolumeSamplingInfo
}

export interface SamplingInfo {
    // Info about "downsampling dimension"
    spatial_downsampling_levels: number[]
    boxes: { [downsampling: number]: SamplingBox }
    time_transformations: TimeTransformation[]
    // e.g. (0, 1, 2) as standard
    original_axis_order: Vector3
    source_axes_units: { [ axis: string ]: string}
}

export interface VolumeSamplingInfo extends SamplingInfo {
    // resolution -> time -> channel_id
    descriptive_statistics: {[resolution: number]: {
        [time: number]: {
            [channel_id: string]: VolumeDescriptiveStatistics}}}
}

export interface VolumeDescriptiveStatistics {
    mean: number
    min: number
    max: number
    std: number
}

export interface TimeTransformation {
    // # to which downsampling level it is applied: can be to specific level, can be to all lvls
    downsampling_level: string | number
    factor: number
}

export interface TimeInfo {
    // #just one kind - range
    kind: string
    start: number
    end: number
    units: string
}

export interface SamplingBox {
    origin: Vector3
    voxel_size: Vector3
    grid_dimensions: Vector3
}

export interface SegmentationLattices {
    segmentation_lattice_ids: number[],
    segmentation_downsamplings: { [lattice: number]: number[] },
}

export interface MeshesMetadata {
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

export interface SegmentationLatticesMetadata {
    segmentation_lattice_ids: number[]
    // #maps lattice id to Sampling Info
    segmentation_sampling_info: { [lattice_id: number]: SamplingInfo }
    // #maps lattice id to channel_ids for that lattice
    channel_ids: { [lattice_id: number]: number[] }
    // #maps lattice id to TimeInfo
    time_info: { [lattice_id: number]: TimeInfo }
}

export interface EntryId {
    source_db_name: string
    source_db_id: string
}

export interface AnnotationMetadata {
    entry_id: EntryId
    details: string | undefined
    name: string
    non_segment_annotation: BiologicalAnnotation | undefined
    volume_channels_annotations: ChannelAnnotation[] | undefined
    segmentation_lattices: SegmentationLatticeInfo[]
}

export interface ChannelAnnotation {
    channel_id: number
    // # with transparency
    color: Vector4
    label: string | undefined
}

export interface Segment {
    id: number
    color: Vector4
    biological_annotation: BiologicalAnnotation
    extra_annotations: BiologicalAnnotation | undefined
}

export interface SegmentationLatticeInfo {
    lattice_id: number
    segment_list: Segment[]
}

export interface BiologicalAnnotation {
    name: string
    external_references: ExternalReference[]
    is_hidden: boolean | undefined
}

export interface ExternalReference {
    id: number, resource: string, accession: string, label: string,
    description: string
}

type Vector3 = [number, number, number];
type Vector4 = [number, number, number, number];