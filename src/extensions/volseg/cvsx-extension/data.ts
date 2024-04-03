
import { Source } from '../common';
import { AnnotationMetadata, GridMetadata, ShapePrimitiveData } from '../new-volumes-and-segmentations/volseg-api/data';

export interface CVSXFilesData {
    volumes?: CVSXVolumeData[]
    latticeSegmentations?: CVSXLatticeSegmentationData[],
    geometricSegmentations?: CVSXGeometricSegmentationData[],
    meshSegmentations?: CVSXMeshSegmentationData[],
    annotation?: AnnotationMetadata,
    metadata?: GridMetadata,
    query: QueryArgs
};

export interface CVSXGeometricSegmentationData {
    segmentationId: string
    timeframeIndex: number
    data: ShapePrimitiveData
}

export interface CVSXMeshSegmentationData {
    segmentationId: string
    timeframeIndex: number
    // segment id, segment data
    data: [string, Uint8Array][]
}

export interface CVSXVolumeData {
    channelId: string
    timeframeIndex: number
    data: Uint8Array
}

export interface CVSXLatticeSegmentationData {
    segmentationId: string
    timeframeIndex: number
    data: Uint8Array
}

export interface QueryArgs {
    entry_id: string,
    source_db: Source,
    time?: number,
    channel_id?: string,
    segmentation_id?: string,
    detail_lvl?: number,
    max_points?: number
}