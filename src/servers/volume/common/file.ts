/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs'
import * as path from 'path'
import * as DataFormat from './data-format'
import { FileHandle } from 'mol-io/common/file-handle';
import { SimpleBuffer } from 'mol-io/common/simple-buffer';

export const IsNativeEndianLittle = new Uint16Array(new Uint8Array([0x12, 0x34]).buffer)[0] === 0x3412;

export async function openRead(filename: string) {
    return new Promise<number>((res, rej) => {
        fs.open(filename, 'r', async (err, file) => {
            if (err) {
                rej(err);
                return;
            }

            try {
                res(file);
            } catch (e) {
                fs.closeSync(file);
            }
        })
    });
}

function makeDir(path: string, root?: string): boolean {
    let dirs = path.split(/\/|\\/g),
        dir = dirs.shift();

    root = (root || '') + dir + '/';

    try { fs.mkdirSync(root); }
    catch (e) {
        if (!fs.statSync(root).isDirectory()) throw new Error(e);
    }

    return !dirs.length || makeDir(dirs.join('/'), root);
}

export function exists(filename: string) {
    return fs.existsSync(filename);
}

export function createFile(filename: string) {
    return new Promise<number>((res, rej) => {
        if (fs.existsSync(filename)) fs.unlinkSync(filename);
        makeDir(path.dirname(filename));
        fs.open(filename, 'w', (err, file) => {
            if (err) rej(err);
            else res(file);
        })
    });
}

const smallBuffer = SimpleBuffer.fromBuffer(new Buffer(8));
export async function writeInt(file: FileHandle, value: number, position: number) {
    smallBuffer.writeInt32LE(value, 0);
    await file.writeBuffer(position, smallBuffer, 4);
}

export interface TypedArrayBufferContext {
    type: DataFormat.ValueType,
    elementByteSize: number,
    readBuffer: Buffer,
    valuesBuffer: Uint8Array,
    values: DataFormat.ValueArray
}

function getElementByteSize(type: DataFormat.ValueType) {
    if (type === DataFormat.ValueType.Float32) return 4;
    if (type === DataFormat.ValueType.Int16) return 2;
    return 1;
}

function makeTypedArray(type: DataFormat.ValueType, buffer: ArrayBuffer): DataFormat.ValueArray {
    if (type === DataFormat.ValueType.Float32) return new Float32Array(buffer);
    if (type === DataFormat.ValueType.Int16) return new Int16Array(buffer);
    return new Int8Array(buffer);
}

export function createTypedArrayBufferContext(size: number, type: DataFormat.ValueType): TypedArrayBufferContext {
    let elementByteSize = getElementByteSize(type);
    let arrayBuffer = new ArrayBuffer(elementByteSize * size);
    let readBuffer = new Buffer(arrayBuffer);
    let valuesBuffer = IsNativeEndianLittle ? arrayBuffer : new ArrayBuffer(elementByteSize * size);
    return {
        type,
        elementByteSize,
        readBuffer,
        valuesBuffer: new Uint8Array(valuesBuffer),
        values: makeTypedArray(type, valuesBuffer)
    };
}

function flipByteOrder(source: SimpleBuffer, target: Uint8Array, byteCount: number, elementByteSize: number, offset: number) {
    for (let i = 0, n = byteCount; i < n; i += elementByteSize) {
        for (let j = 0; j < elementByteSize; j++) {
            target[offset + i + elementByteSize - j - 1] = source[offset + i + j];
        }
    }
}

export async function readTypedArray(ctx: TypedArrayBufferContext, file: FileHandle, position: number, count: number, valueOffset: number, littleEndian?: boolean) {
    let byteCount = ctx.elementByteSize * count;
    let byteOffset = ctx.elementByteSize * valueOffset;

    await file.readBuffer(position, ctx.readBuffer, byteCount, byteOffset);
    if (ctx.elementByteSize > 1 && ((littleEndian !== void 0 && littleEndian !== IsNativeEndianLittle) || !IsNativeEndianLittle)) {
        // fix the endian
        flipByteOrder(ctx.readBuffer, ctx.valuesBuffer, byteCount, ctx.elementByteSize, byteOffset);
    }
    return ctx.values;
}

export function ensureLittleEndian(source: SimpleBuffer, target: SimpleBuffer, byteCount: number, elementByteSize: number, offset: number) {
    if (IsNativeEndianLittle) return;
    if (!byteCount || elementByteSize <= 1) return;
    flipByteOrder(source, target, byteCount, elementByteSize, offset);
}