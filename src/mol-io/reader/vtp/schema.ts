/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

export interface VtpDataArrayDescriptor {
    name: string;
    /** VTK type string: "Float32", "Float64", "Int32", "Int64", "UInt32", etc. */
    type: string;
    numberOfComponents: number;
    /** Byte offset into the appended data section (after the `_` marker) */
    offset: number;
    rangeMin: number;
    rangeMax: number;
}

export interface VtpScalarArray {
    desc: VtpDataArrayDescriptor;
    /** Scalar values (one per point or cell), always stored as Float64 */
    values: Float64Array;
}

export interface VtpFile {
    numberOfPoints: number;
    /** Number of cells before fan-triangulation */
    numberOfCells: number;
    /** Vertex positions: [x0,y0,z0, x1,y1,z1, ...], length = 3 * numberOfPoints */
    positions: Float32Array;
    /**
     * Triangle connectivity after fan-triangulation: [v0,v1,v2, ...],
     * length = 3 * numberOfTriangles.
     */
    connectivity: Int32Array;
    /** Total number of triangles (≥ numberOfCells for quads). */
    numberOfTriangles: number;
    /** Per-vertex (PointData) scalar arrays. Key = array name. */
    pointData: Map<string, VtpScalarArray>;
    /** Per-cell (CellData) scalar arrays. Key = array name. */
    cellData: Map<string, VtpScalarArray>;
}
