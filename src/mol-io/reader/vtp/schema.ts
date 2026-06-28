/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { Column } from '../../../mol-data/db';

/** Scalar element types supported by VTK DataArray elements. */
export type VtkDataType =
    | 'Float32' | 'Float64'
    | 'Int8' | 'UInt8'
    | 'Int16' | 'UInt16'
    | 'Int32' | 'UInt32'
    | 'Int64' | 'UInt64';

const VtkDataTypeSet: ReadonlySet<string> = new Set<VtkDataType>([
    'Float32', 'Float64', 'Int8', 'UInt8', 'Int16', 'UInt16', 'Int32', 'UInt32', 'Int64', 'UInt64',
]);

/** Narrow a raw VTK `type` attribute to a VtkDataType, falling back to 'Float32' if unrecognized. */
export function toVtkDataType(type: string | undefined): VtkDataType {
    return type !== undefined && VtkDataTypeSet.has(type) ? type as VtkDataType : 'Float32';
}

/**
 * Decoded data array surfaced to consumers. Only the metadata downstream code needs is retained;
 * the transient parsing descriptor (raw base64/ascii source, byte offsets, encoding format) is kept
 * private to the parser so its potentially large source strings can be garbage-collected after decode.
 */
export interface VtpScalarArray {
    readonly name: string;
    readonly type: VtkDataType;
    readonly numberOfComponents: number;
    /** Value range as declared in the XML (RangeMin/RangeMax). undefined when the writer omitted it (many do). */
    readonly rangeMin: number | undefined;
    readonly rangeMax: number | undefined;
    /**
     * Flat decoded values as a Column. For scalar arrays (numberOfComponents=1) rowCount = nElems.
     * For multi-component arrays rowCount = nElems * numberOfComponents.
     * Backed by a typed array for now; the Column interface lets us defer full decoding later
     * without changing this contract. Int64/UInt64 are represented as Float64-backed values.
     */
    readonly values: Column<number>;
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
