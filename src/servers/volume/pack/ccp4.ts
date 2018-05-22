/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as File from '../common/file'
import * as DataFormat from '../common/data-format'

export const enum Mode { Int8 = 0, Int16 = 1, Float32 = 2 }

export interface Header {
    name: string,
    mode: Mode,
    grid: number[], // grid is converted to the axis order!!
    axisOrder: number[],
    extent: number[],
    origin: number[],
    spacegroupNumber: number,
    cellSize: number[],
    cellAngles: number[],
    littleEndian: boolean,
    dataOffset: number
}

/** Represents a circular buffer for 2 * blockSize layers */
export interface SliceBuffer {
    buffer: File.TypedArrayBufferContext,
    sliceCapacity: number,
    slicesRead: number,

    values: DataFormat.ValueArray,
    sliceCount: number,

    /** Have all the input slice been read? */
    isFinished: boolean
}

export interface Data {
    header: Header,
    file: number,
    slices: SliceBuffer
}

export function getValueType(header: Header) {
    if (header.mode === Mode.Float32) return DataFormat.ValueType.Float32;
    if (header.mode === Mode.Int16) return DataFormat.ValueType.Int16;
    return DataFormat.ValueType.Int8;
}

export function assignSliceBuffer(data: Data, blockSize: number) {
    const { extent } = data.header;
    const valueType = getValueType(data.header);
    const sliceSize = extent[0] * extent[1] * DataFormat.getValueByteSize(valueType);
    const sliceCapacity = Math.max(1, Math.floor(Math.min(64 * 1024 * 1024, sliceSize * extent[2]) / sliceSize));
    const buffer = File.createTypedArrayBufferContext(sliceCapacity * extent[0] * extent[1], valueType);
    data.slices = {
        buffer,
        sliceCapacity,
        slicesRead: 0,
        values: buffer.values,
        sliceCount: 0,
        isFinished: false
    };
}

function compareProp(a: any, b: any) {
    if (a instanceof Array && b instanceof Array) {
        if (a.length !== b.length) return false;
        for (let i = 0; i < a.length; i++) {
            if (a[i] !== b[i]) return false;
        }
        return true;
    }
    return a === b;
}

export function compareHeaders(a: Header, b: Header) {
    for (const p of ['grid', 'axisOrder', 'extent', 'origin', 'spacegroupNumber', 'cellSize', 'cellAngles', 'mode']) {
        if (!compareProp((a as any)[p], (b as any)[p])) return false;
    }
    return true;
}

function getArray(r: (offset: number) => number, offset: number, count: number) {
    const ret: number[] = [];
    for (let i = 0; i < count; i++) {
        ret[i] = r(offset + i);
    }
    return ret;
}

async function readHeader(name: string, file: number) {
    const headerSize = 1024;
    const { buffer: data } = await File.readBuffer(file, 0, headerSize);

    let littleEndian = true;

    let mode = data.readInt32LE(3 * 4);
    if (mode < 0 || mode > 2) {
        littleEndian = false;
        mode = data.readInt32BE(3 * 4, true);
        if (mode < 0 || mode > 2) {
            throw Error('Only CCP4 modes 0, 1, and 2 are supported.');
        }
    }

    const readInt = littleEndian ? (o: number) => data.readInt32LE(o * 4) : (o: number) => data.readInt32BE(o * 4);
    const readFloat = littleEndian ? (o: number) => data.readFloatLE(o * 4) : (o: number) => data.readFloatBE(o * 4);

    const origin2k = getArray(readFloat, 49, 3);
    const nxyzStart = getArray(readInt, 4, 3);
    const header: Header = {
        name,
        mode,
        grid: getArray(readInt, 7, 3),
        axisOrder: getArray(readInt, 16, 3).map(i => i - 1),
        extent: getArray(readInt, 0, 3),
        origin: origin2k[0] === 0.0 && origin2k[1] === 0.0 && origin2k[2] === 0.0 ? nxyzStart : origin2k,
        spacegroupNumber: readInt(22),
        cellSize: getArray(readFloat, 10, 3),
        cellAngles: getArray(readFloat, 13, 3),
        // mean: readFloat(21),
        littleEndian,
        dataOffset: headerSize + readInt(23) /* symBytes */
    };
    // "normalize" the grid axis order
    header.grid = [header.grid[header.axisOrder[0]], header.grid[header.axisOrder[1]], header.grid[header.axisOrder[2]]];
    return header;
}

export async function readSlices(data: Data) {
    const { slices, header } = data;
    if (slices.isFinished) {
        return;
    }

    const { extent } = header;
    const sliceSize = extent[0] * extent[1];
    const sliceByteOffset = slices.buffer.elementByteSize * sliceSize * slices.slicesRead;
    const sliceCount = Math.min(slices.sliceCapacity, extent[2] - slices.slicesRead);
    const sliceByteCount = sliceCount * sliceSize;

    await File.readTypedArray(slices.buffer, data.file, header.dataOffset + sliceByteOffset, sliceByteCount, 0, header.littleEndian);
    slices.slicesRead += sliceCount;
    slices.sliceCount = sliceCount;

    if (slices.slicesRead >= extent[2]) {
        slices.isFinished = true;
    }
}

export async function open(name: string, filename: string): Promise<Data> {
    const file = await File.openRead(filename);
    const header = await readHeader(name, file);
    return {
        header,
        file,
        slices: void 0 as any
    };
}