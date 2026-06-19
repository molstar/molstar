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
    /** Encoding format of this DataArray element. */
    readonly format: 'ascii' | 'binary' | 'appended';
    /** Byte/char offset into the appended data section (after the `_` marker). -1 for inline/ascii. */
    readonly offset: number;
    readonly rangeMin: number;
    readonly rangeMax: number;
    /** Base64 content of the DataArray element (inline binary format only). */
    readonly inlineBase64?: string;
    /** Raw text content of the DataArray element (ascii format only). */
    readonly asciiText?: string;
}

export interface VtpScalarArray {
    readonly desc: VtpDataArrayDescriptor;
    /**
     * Flat decoded values. For scalar arrays (numberOfComponents=1) length = nElems.
     * For multi-component arrays length = nElems * numberOfComponents.
     */
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
    /** Per-vertex (PointData) arrays. Key = array name. */
    readonly pointData: ReadonlyMap<string, VtpScalarArray>;
    /** Per-cell (CellData) arrays. Key = array name. */
    readonly cellData: ReadonlyMap<string, VtpScalarArray>;
}
