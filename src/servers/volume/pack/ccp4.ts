/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as File from '../common/file'
import * as DataFormat from '../common/data-format'
import { FileHandle } from 'mol-io/common/file-handle';
import { readCcp4Header } from 'mol-io/reader/ccp4/parser';

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
    file: FileHandle,
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

async function readHeader(name: string, file: FileHandle) {
    const { header: ccp4Header, littleEndian } = await readCcp4Header(file)

    const origin2k = [ccp4Header.originX, ccp4Header.originY, ccp4Header.originZ];
    const nxyzStart = [ccp4Header.NCSTART, ccp4Header.NRSTART, ccp4Header.NSSTART];
    const header: Header = {
        name,
        mode: ccp4Header.MODE,
        grid: [ccp4Header.NX, ccp4Header.NY, ccp4Header.NZ],
        axisOrder: [ccp4Header.MAPC, ccp4Header.MAPR, ccp4Header.MAPS].map(i => i - 1),
        extent: [ccp4Header.NC, ccp4Header.NR, ccp4Header.NS],
        origin: origin2k[0] === 0.0 && origin2k[1] === 0.0 && origin2k[2] === 0.0 ? nxyzStart : origin2k,
        spacegroupNumber: ccp4Header.ISPG,
        cellSize: [ccp4Header.xLength, ccp4Header.yLength, ccp4Header.zLength],
        cellAngles: [ccp4Header.alpha, ccp4Header.beta, ccp4Header.gamma],
        // mean: readFloat(21),
        littleEndian,
        dataOffset: 256 * 4 + ccp4Header.NSYMBT /* symBytes */
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
    const descriptor = await File.openRead(filename);
    const file = FileHandle.fromDescriptor(descriptor)
    const header = await readHeader(name, file);
    return {
        header,
        file,
        slices: void 0 as any
    };
}