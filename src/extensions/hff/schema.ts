/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Typed representation of the EMDB-SFF (Segmentation File Format) data model
 * as serialized to HDF5 (.hff / .h5 / .hdf5). Schema reference:
 *   https://emdb-empiar.github.io/EMDB-SFF/
 *   https://github.com/emdb-empiar/sfftk-rw  (writer authoritative)
 *
 * v1 scope: mesh_list segments. Lattices and shape primitives are typed
 * but not parsed yet — those readers will be added later.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

export type SffMode =
    | 'int8' | 'uint8'
    | 'int16' | 'uint16'
    | 'int32' | 'uint32'
    | 'int64' | 'uint64'
    | 'float32' | 'float64';

export type SffEndianness = 'little' | 'big';

/** RGBA in [0..1] floats. */
export type SffColour = [number, number, number, number];

export interface SffEncodedSequence {
    /** Number of (x,y,z) triples. */
    count: number;
    mode: SffMode;
    endianness: SffEndianness;
    /** Decoded (n,3) flat array; native endianness, length = count*3. */
    data: Float32Array | Float64Array | Int32Array | Uint32Array | Int16Array | Uint16Array | Int8Array | Uint8Array;
}

export interface SffMesh {
    id: number;
    vertices: SffEncodedSequence;
    triangles: SffEncodedSequence;
    normals?: SffEncodedSequence;
    /** Index into SffData.transforms to apply to vertices/normals (per spec). */
    transformId?: number;
}

export interface SffBiologicalAnnotation {
    name?: string;
    description?: string;
    numberOfInstances?: number;
}

export interface SffSegment {
    id: number;
    parentId: number;
    colour: SffColour;
    biologicalAnnotation?: SffBiologicalAnnotation;
    meshes: SffMesh[];
    /** Lattice/three_d_volume/shape data are read into raw blobs for now. */
    threeDVolume?: SffThreeDVolume;
    shapePrimitives?: unknown;
}

export interface SffThreeDVolume {
    latticeId: number;
    value: number;
    transformId?: number;
}

export interface SffTransform {
    id: number;
    rows: number;
    cols: number;
    /** Row-major matrix data (rows*cols floats). */
    data: Float64Array;
}

export type SffPrimaryDescriptor = 'three_d_volume' | 'mesh_list' | 'shape_primitive_list';

export interface SffData {
    version?: string;
    name?: string;
    details?: string;
    primaryDescriptor?: SffPrimaryDescriptor;
    transforms: SffTransform[];
    segments: SffSegment[];
}
