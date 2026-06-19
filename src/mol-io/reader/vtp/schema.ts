/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

export interface VtpDataArrayDescriptor {
    readonly name: string;
    /** VTK type string: "Float32", "Float64", "Int32", "Int64", "UInt32", etc. */
    readonly type: string;
    readonly numberOfComponents: number;
    /** Byte/char offset into the appended data section (after the `_` marker). -1 for inline binary. */
    readonly offset: number;
    readonly rangeMin: number;
    readonly rangeMax: number;
    /** Base64 content of the DataArray element (inline binary format only). */
    readonly inlineBase64?: string;
}

export interface VtpScalarArray {
    readonly desc: VtpDataArrayDescriptor;
    readonly values: Float64Array;
}

export interface VtpFile {
    readonly numberOfPoints: number;
    readonly numberOfCells: number;
    /** Vertex positions: [x0,y0,z0, x1,y1,z1, ...], length = 3 * numberOfPoints */
    readonly positions: Float32Array;
    /** Triangle connectivity after fan-triangulation: [v0,v1,v2, ...], length = 3 * numberOfTriangles */
    readonly connectivity: Int32Array;
    readonly numberOfTriangles: number;
    /** Maps triangle index → originating VTK cell index, length = numberOfTriangles */
    readonly triangleCellIndex: Int32Array;
    /** Per-vertex (PointData) scalar arrays. Key = array name. */
    readonly pointData: ReadonlyMap<string, VtpScalarArray>;
    /** Per-cell (CellData) scalar arrays. Key = array name. */
    readonly cellData: ReadonlyMap<string, VtpScalarArray>;
}
