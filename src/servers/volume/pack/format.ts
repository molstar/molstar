/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as File from '../common/file';
import { FileHandle } from '../../../mol-io/common/file-handle';
import { Ccp4Provider } from './format/ccp4';
import { TypedArrayBufferContext, TypedArrayValueArray, TypedArrayValueType, getElementByteSize, createTypedArrayBufferContext } from '../../../mol-io/common/typed-array';
import { Dsn6Provider } from './format/dsn6';

export interface Header {
    name: string,
    valueType: TypedArrayValueType,
    grid: number[], // grid is converted to the axis order!!
    axisOrder: number[],
    extent: number[],
    origin: number[],
    spacegroupNumber: number,
    cellSize: number[],
    cellAngles: number[],
    littleEndian: boolean,
    dataOffset: number
    originalHeader: unknown // TODO
}

/** Represents a circular buffer for 2 * blockSize layers */
export interface SliceBuffer {
    buffer: TypedArrayBufferContext,
    maxBlockBytes: number
    sliceCapacity: number,
    slicesRead: number,

    values: TypedArrayValueArray,
    sliceCount: number,

    /** Have all the input slice been read? */
    isFinished: boolean
}

export interface Data {
    header: Header,
    file: FileHandle,
    slices: SliceBuffer
}

export interface Provider {
    readHeader: (name: string, file: FileHandle) => Promise<Header>,
    readSlices: (data: Data) => Promise<void>
}

export interface Context {
    data: Data,
    provider: Provider
}

export function assignSliceBuffer(data: Data, blockSizeInMB: number) {
    const { extent, valueType } = data.header;
    const maxBlockBytes = blockSizeInMB * 1024 * 1024;
    const sliceSize = extent[0] * extent[1] * getElementByteSize(valueType);
    const sliceCapacity = Math.max(1, Math.floor(Math.min(maxBlockBytes, sliceSize * extent[2]) / sliceSize));
    const buffer = createTypedArrayBufferContext(sliceCapacity * extent[0] * extent[1], valueType);
    data.slices = {
        buffer,
        maxBlockBytes,
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

export type Type = 'ccp4' | 'dsn6'

export function getProviderFromType(type: Type): Provider {
    switch (type) {
        case 'ccp4': return Ccp4Provider;
        case 'dsn6': return Dsn6Provider;
    }
}

export async function open(name: string, filename: string, type: Type): Promise<Context> {
    const provider = getProviderFromType(type);
    const descriptor = await File.openRead(filename);
    const file = FileHandle.fromDescriptor(descriptor, filename);
    const header = await provider.readHeader(name, file);
    const data = { header, file, slices: void 0 as any };
    return { data, provider };
}