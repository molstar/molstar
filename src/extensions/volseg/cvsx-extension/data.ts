import { AnnotationMetadata, GridMetadata, ShapePrimitiveData } from ../ new- volumes - and - segmentations / volseg - api / data

export interface CVSXFilesData {
    volumes?: CVSXVolumeData[]
    latticeSegmentations?: CVSXLatticeSegmentationData[],
    // TODO: other
    geometricSegmentations?: ShapePrimitiveData,
    meshSegmentations?: [string, Uint8Array][],
    annotation?: AnnotationMetadata,
    metadata?: GridMetadata,
    query: QueryArgs
};

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

// export interface QueryMetadata {
//     subquery_types: string[]
//     args: QueryArgs
// }

interface QueryArgs {
    entry_id: string,
    source_db: string,
    time?: number,
    channel_id?: string,
    segmentation_id?: string,
    detail_lvl?: number,
    max_points?: number
}