import { AnnotationMetadata, GridMetadata, ShapePrimitiveData } from ../ new- volumes - and - segmentations / volseg - api / data

export interface CVSXFilesData {
    // parsed everything
    volume?: Uint8Array
    latticeSegmentation?: Uint8Array,
    geometricSegmentation?: ShapePrimitiveData,
    meshSegmentation?: [string, Uint8Array][],
    annotation?: AnnotationMetadata,
    metadata?: GridMetadata,
    query: QueryMetadata
};

export interface QueryMetadata {
    subquery_types: string[]
    args: QueryArgs
}

interface QueryArgs {
    entry_id: string,
    source_db: string,
    time?: number,
    channel_id?: string,
    segmentation_id?: string,
    detail_lvl?: number,
    max_points?: number
}