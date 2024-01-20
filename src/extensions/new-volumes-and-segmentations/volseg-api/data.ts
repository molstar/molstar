/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

export interface Metadata {
    grid: GridMetadata,
    annotation?: AnnotationMetadata
}

export interface GridMetadata {
    entry_id: EntryId
    volumes: VolumesMetadata
    segmentation_lattices?: SegmentationLatticesMetadata
    segmentation_meshes?: MeshSegmentationSetsMetadata
    geometric_segmentation?: GeometricSegmentationSetsMetadata
}

export interface MeshSegmentationSetsMetadata {
    segmentation_ids: string[]
    segmentation_metadata: { [segmentation_id: string]: MeshesMetadata }
    // maps set ids to time info
    time_info: { [segmentation_id: string]: TimeInfo }
}

export interface MeshesMetadata {
    segmentation_mesh_set_id: string
    // maps timeframe index to MeshesTimeframeMetadata with mesh comp num
    mesh_timeframes: { [timeframe_index: number]: MeshComponentNumbers }
    detail_lvl_to_fraction: {
        [lvl: number]: number
    }
}

export interface MeshComponentNumbers {
    segment_ids?: {
        [segId: number]: {
            detail_lvls: {
                [detail: number]: {
                    mesh_ids: {
                        [meshId: number]: {
                            num_triangles: number
                            num_vertices: number
                            num_normals?: number
                        }
                    }
                }
            }
        }
    }
}

export interface EntryId {
    source_db_name: string
    source_db_id: string
}

export interface GeometricSegmentationSetsMetadata {
    segmentation_ids: string[]
    // maps set ids to time info
    time_info: { [segmentation_id: string]: TimeInfo }
}

export interface VolumesMetadata {
    channel_ids: string[]
    // Values of time dimension
    time_info: TimeInfo
    volume_sampling_info: VolumeSamplingInfo
}

export interface SamplingInfo {
    // Info about 'downsampling dimension'
    spatial_downsampling_levels: number[]
    boxes: { [downsampling: number]: SamplingBox }
    time_transformations: TimeTransformation[]
    // e.g. (0, 1, 2) as standard
    original_axis_order: Vector3
    source_axes_units: { [ axis: string ]: string}
}

export interface VolumeSamplingInfo extends SamplingInfo {
    // resolution -> time -> channel_id
    descriptive_statistics: {
        [resolution: number]: {
            [time: number]: {
                [channel_id: string]: VolumeDescriptiveStatistics
            }
        }
    }
}

export interface VolumeDescriptiveStatistics {
    mean: number
    min: number
    max: number
    std: number
}

export interface TimeTransformation {
    // # to which downsampling level it is applied: can be to specific level, can be to all lvls
    downsampling_level: 'all' | number
    factor: number
}

export interface TimeInfo {
    // just one kind - range
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

export interface SegmentationLatticesMetadata {
    // e.g. label groups (Cell, Chromosomes)
    segmentation_lattice_ids: string[]
    segmentation_sampling_info: { [lattice_id: string]: SamplingInfo }
    // #maps lattice id to TimeInfo
    time_info: { [segmentation_id: string]: TimeInfo }
}

export interface EntryId {
    source_db_name: string
    source_db_id: string
}

export interface AnnotationMetadata {
    name?: string
    entry_id: EntryId
    // id => DescriptionData
    descriptions: { [id: string]: DescriptionData }
    // NOTE: on frontend, segment key = `${kind}:{segmentation_id}:{segment_id}`
    annotations: SegmentAnnotationData[]
    details?: string
    volume_channels_annotations?: ChannelAnnotation[]
}

export interface ChannelAnnotation {
    channel_id: string
    color: Vector4
    label?: string
}

export interface DescriptionData {
    id: string
    target_kind: 'lattice' | 'mesh' | 'primitive' | 'entry'
    target_id?: TargetId
    name?: string
    external_references?: ExternalReference[]
    is_hidden?: boolean
    time?: number | number[] | Vector2[]

    description?: Description
    metadata?: { [key: string]: any}
}

export interface TargetId {
    segmentation_id: string
    segment_id: number
}

export interface Description {
    description_format: 'text' | 'markdown'
    description_text: string
}


export interface SegmentAnnotationData {
    id: string
    segment_kind: 'lattice' | 'mesh' | 'primitive'
    segment_id: number
    segmentation_id: string
    color?: Vector4
    time?: number | number[] | Vector2[]
}

export interface ExternalReference {
    id: number
    resource?: string
    accession?: string
    label?: string,
    description?: string
}

// export interface Segment {
//     id: number
//     color: Vector4
//     biological_annotation: BiologicalAnnotation
//     extra_annotations: BiologicalAnnotation | undefined
// }

// export interface SegmentationLatticeInfo {
//     lattice_id: number
//     segment_list: Segment[]
// }

// export interface BiologicalAnnotation {
//     name: string
//     external_references: ExternalReference[]
//     is_hidden: boolean | undefined
// }

// export interface ExternalReference {
//     id: number, resource: string, accession: string, label: string,
//     description: string
// }

type Vector2 = [number, number];
type Vector3 = [number, number, number];
type Vector4 = [number, number, number, number];