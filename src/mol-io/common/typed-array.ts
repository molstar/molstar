/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { FileHandle } from '../../mol-io/common/file-handle';
import { SimpleBuffer } from '../../mol-io/common/simple-buffer';

export type TypedArrayValueType = 'float32' | 'int8' | 'int16' | 'uint16'

export namespace TypedArrayValueType {
    export const Float32: TypedArrayValueType = 'float32';
    export const Int8: TypedArrayValueType = 'int8';
    export const Int16: TypedArrayValueType = 'int16';
    export const Uint16: TypedArrayValueType = 'uint16';
}

export type TypedArrayValueArray = Float32Array | Int8Array | Int16Array | Uint16Array

export interface TypedArrayBufferContext {
    type: TypedArrayValueType,
    elementByteSize: number,
    readBuffer: SimpleBuffer,
    valuesBuffer: Uint8Array,
    values: TypedArrayValueArray
}

export function getElementByteSize(type: TypedArrayValueType) {
    if (type === TypedArrayValueType.Float32) return 4;
    if (type === TypedArrayValueType.Int16) return 2;
    if (type === TypedArrayValueType.Uint16) return 2;
    return 1;
}

export function makeTypedArray(type: TypedArrayValueType, buffer: ArrayBuffer, byteOffset = 0, length?: number): TypedArrayValueArray {
    if (type === TypedArrayValueType.Float32) return new Float32Array(buffer, byteOffset, length);
    if (type === TypedArrayValueType.Int16) return new Int16Array(buffer, byteOffset, length);
    if (type === TypedArrayValueType.Uint16) return new Uint16Array(buffer, byteOffset, length);
    return new Int8Array(buffer, byteOffset, length);
}

export function createTypedArray(type: TypedArrayValueType, size: number) {
    switch (type) {
        case TypedArrayValueType.Float32: return new Float32Array(new ArrayBuffer(4 * size));
        case TypedArrayValueType.Int8: return new Int8Array(new ArrayBuffer(1 * size));
        case TypedArrayValueType.Int16: return new Int16Array(new ArrayBuffer(2 * size));
        case TypedArrayValueType.Uint16: return new Uint16Array(new ArrayBuffer(2 * size));
    }
    throw Error(`${type} is not a supported value format.`);
}

export function createTypedArrayBufferContext(size: number, type: TypedArrayValueType): TypedArrayBufferContext {
    let elementByteSize = getElementByteSize(type);
    let arrayBuffer = new ArrayBuffer(elementByteSize * size);
    let readBuffer = SimpleBuffer.fromArrayBuffer(arrayBuffer);
    let valuesBuffer = SimpleBuffer.IsNativeEndianLittle ? arrayBuffer : new ArrayBuffer(elementByteSize * size);
    return {
        type,
        elementByteSize,
        readBuffer,
        valuesBuffer: new Uint8Array(valuesBuffer),
        values: makeTypedArray(type, valuesBuffer)
    };
}

export async function readTypedArray(ctx: TypedArrayBufferContext, file: FileHandle, position: number, byteCount: number, valueByteOffset: number, littleEndian?: boolean) {
    await file.readBuffer(position, ctx.readBuffer, byteCount, valueByteOffset);
    if (ctx.elementByteSize > 1 && ((littleEndian !== void 0 && littleEndian !== SimpleBuffer.IsNativeEndianLittle) || !SimpleBuffer.IsNativeEndianLittle)) {
        // fix the endian
        SimpleBuffer.flipByteOrder(ctx.readBuffer, ctx.valuesBuffer, byteCount, ctx.elementByteSize, valueByteOffset);
    }
    return ctx.values;
}
