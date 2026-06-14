/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '../../../mol-task';
import { ReaderResult as Result } from '../result';

// EM file format specification:
// https://github.com/uermel/Artiatomi/wiki/EM-(file-format)

/** EM file header size in bytes. */
export const EmHeader_Size = 512;

/** Machine coding values stored in byte 0 of the EM header. */
export const EmMachineCoding = {
    OS9: 0,
    VAX: 1,
    Convex: 2,
    SGI: 3,
    Mac: 5,
    PC: 6,
} as const;

/** Data type coding values stored in byte 3 of the EM header. */
export const EmDataType = {
    Byte: 1, // 1 byte  (signed int8)
    Short: 2, // 2 bytes (signed int16)
    Int: 4, // 4 bytes (signed int32)
    Float: 5, // 4 bytes (float32)
    Complex: 8, // 8 bytes (two float32 values)
    Double: 9, // 8 bytes (float64)
} as const;

export interface ArtiatomiEmHeader {
    readonly machineCoding: number
    readonly dataType: number
    readonly dimX: number
    readonly dimY: number
    readonly dimZ: number
    readonly isLittleEndian: boolean
}

/** Typed array union for EM data. For Complex, the array has length 2*N with interleaved (real, imaginary) pairs. */
export type EmData = Int8Array | Int16Array | Int32Array | Float32Array | Float64Array;

export interface ArtiatomiEmFile {
    readonly header: ArtiatomiEmHeader
    /**
     * Raw numeric data stored in the appropriate typed array for the data type:
     * - Byte    → Int8Array,   length = N
     * - Short   → Int16Array,  length = N
     * - Int     → Int32Array,  length = N
     * - Float   → Float32Array, length = N
     * - Complex → Float32Array, length = 2*N, interleaved as (r0,i0, r1,i1, ...)
     * - Double  → Float64Array, length = N
     * where N = dimX * dimY * dimZ and X is the fastest-varying dimension.
     */
    readonly data: EmData
}

function isLittleEndian(machineCoding: number): boolean {
    // VAX (1) and PC (6) write little-endian data; all other machine codes are big-endian.
    return machineCoding === EmMachineCoding.VAX || machineCoding === EmMachineCoding.PC;
}

function bytesPerElement(dataType: number): number {
    switch (dataType) {
        case EmDataType.Byte: return 1;
        case EmDataType.Short: return 2;
        case EmDataType.Int: return 4;
        case EmDataType.Float: return 4;
        case EmDataType.Complex: return 8;
        case EmDataType.Double: return 8;
        default: return -1;
    }
}

async function parseInternal(data: Uint8Array): Promise<Result<ArtiatomiEmFile>> {
    if (data.length < EmHeader_Size) {
        return Result.error<ArtiatomiEmFile>(`EM file too short to contain a valid header (${data.length} < ${EmHeader_Size}).`);
    }

    const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
    const machineCoding = view.getUint8(0);
    const dataType = view.getUint8(3);
    const le = isLittleEndian(machineCoding);

    const dimX = view.getInt32(4, le);
    const dimY = view.getInt32(8, le);
    const dimZ = view.getInt32(12, le);

    if (dimX <= 0 || dimY <= 0 || dimZ <= 0) {
        return Result.error<ArtiatomiEmFile>(`Invalid EM file dimensions: ${dimX} x ${dimY} x ${dimZ}.`);
    }

    const bpe = bytesPerElement(dataType);
    if (bpe < 0) {
        return Result.error<ArtiatomiEmFile>(`Unsupported EM data type: ${dataType}.`);
    }

    const elementCount = dimX * dimY * dimZ;
    const expectedSize = EmHeader_Size + elementCount * bpe;
    if (data.length < expectedSize) {
        return Result.error<ArtiatomiEmFile>(`EM file too short: expected ${expectedSize} bytes, got ${data.length}.`);
    }

    const offset = EmHeader_Size;
    let typedData: EmData;

    switch (dataType) {
        case EmDataType.Byte: {
            const arr = new Int8Array(elementCount);
            for (let i = 0; i < elementCount; i++) arr[i] = view.getInt8(offset + i);
            typedData = arr;
            break;
        }
        case EmDataType.Short: {
            const arr = new Int16Array(elementCount);
            for (let i = 0; i < elementCount; i++) arr[i] = view.getInt16(offset + i * 2, le);
            typedData = arr;
            break;
        }
        case EmDataType.Int: {
            const arr = new Int32Array(elementCount);
            for (let i = 0; i < elementCount; i++) arr[i] = view.getInt32(offset + i * 4, le);
            typedData = arr;
            break;
        }
        case EmDataType.Float: {
            const arr = new Float32Array(elementCount);
            for (let i = 0; i < elementCount; i++) arr[i] = view.getFloat32(offset + i * 4, le);
            typedData = arr;
            break;
        }
        case EmDataType.Complex: {
            // Interleaved (real, imaginary) pairs: length = 2 * elementCount.
            const arr = new Float32Array(elementCount * 2);
            for (let i = 0; i < elementCount; i++) {
                arr[i * 2] = view.getFloat32(offset + i * 8, le); // real
                arr[i * 2 + 1] = view.getFloat32(offset + i * 8 + 4, le); // imaginary
            }
            typedData = arr;
            break;
        }
        case EmDataType.Double: {
            const arr = new Float64Array(elementCount);
            for (let i = 0; i < elementCount; i++) arr[i] = view.getFloat64(offset + i * 8, le);
            typedData = arr;
            break;
        }
        default:
            return Result.error<ArtiatomiEmFile>(`Unsupported EM data type: ${dataType}.`);
    }

    const header: ArtiatomiEmHeader = { machineCoding, dataType, dimX, dimY, dimZ, isLittleEndian: le };
    return Result.success<ArtiatomiEmFile>({ header, data: typedData });
}

export function parseArtiatomiEm(data: Uint8Array) {
    return Task.create<Result<ArtiatomiEmFile>>('Parse Artiatomi EM', async () => {
        return await parseInternal(data);
    });
}
