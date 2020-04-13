/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defaults } from '../../mol-util';

export interface SimpleBuffer extends Uint8Array {
    readInt8: (offset: number) => number
    readUInt8: (offset: number) => number

    writeInt8: (value: number, offset: number) => void
    writeUInt8: (value: number, offset: number) => void

    readInt16LE: (offset: number) => number
    readInt32LE: (offset: number) => number
    readUInt16LE: (offset: number) => number
    readUInt32LE: (offset: number) => number
    readFloatLE: (offset: number) => number
    readDoubleLE: (offset: number) => number

    writeInt16LE: (value: number, offset: number) => void
    writeInt32LE: (value: number, offset: number) => void
    writeUInt16LE: (value: number, offset: number) => void
    writeUInt32LE: (value: number, offset: number) => void
    writeFloatLE: (value: number, offset: number) => void
    writeDoubleLE: (value: number, offset: number) => void

    readInt16BE: (offset: number) => number
    readInt32BE: (offset: number) => number
    readUInt16BE: (offset: number) => number
    readUInt32BE: (offset: number) => number
    readFloatBE: (offset: number) => number
    readDoubleBE: (offset: number) => number

    writeInt16BE: (value: number, offset: number) => void
    writeInt32BE: (value: number, offset: number) => void
    writeUInt16BE: (value: number, offset: number) => void
    writeUInt32BE: (value: number, offset: number) => void
    writeFloatBE: (value: number, offset: number) => void
    writeDoubleBE: (value: number, offset: number) => void

    copy: (targetBuffer: Uint8Array, targetStart?: number, sourceStart?: number, sourceEnd?: number) => number
}

export namespace SimpleBuffer {
    export function fromUint8Array(array: Uint8Array): SimpleBuffer {
        const dv = new DataView(array.buffer);
        return Object.assign(array.subarray(0), {
            readInt8: (offset: number) => dv.getInt8(offset),
            readUInt8: (offset: number) => dv.getUint8(offset),
            writeInt8: (value: number, offset: number) => dv.setInt8(offset, value),
            writeUInt8: (value: number, offset: number) => dv.setUint8(offset, value),

            readInt16LE: (offset: number) => dv.getInt16(offset, true),
            readInt32LE: (offset: number) => dv.getInt32(offset, true),
            readUInt16LE: (offset: number) => dv.getUint16(offset, true),
            readUInt32LE: (offset: number) => dv.getUint32(offset, true),
            readFloatLE: (offset: number) => dv.getFloat32(offset, true),
            readDoubleLE: (offset: number) => dv.getFloat64(offset, true),

            writeInt16LE: (value: number, offset: number) => dv.setInt16(offset, value, true),
            writeInt32LE: (value: number, offset: number) => dv.setInt32(offset, value, true),
            writeUInt16LE: (value: number, offset: number) => dv.setUint16(offset, value, true),
            writeUInt32LE: (value: number, offset: number) => dv.setUint32(offset, value, true),
            writeFloatLE: (value: number, offset: number) => dv.setFloat32(offset, value, true),
            writeDoubleLE: (value: number, offset: number) => dv.setFloat64(offset, value, true),

            readInt16BE: (offset: number) => dv.getInt16(offset, false),
            readInt32BE: (offset: number) => dv.getInt32(offset, false),
            readUInt16BE: (offset: number) => dv.getUint16(offset, false),
            readUInt32BE: (offset: number) => dv.getUint32(offset, false),
            readFloatBE: (offset: number) => dv.getFloat32(offset, false),
            readDoubleBE: (offset: number) => dv.getFloat64(offset, false),

            writeInt16BE: (value: number, offset: number) => dv.setInt16(offset, value, false),
            writeInt32BE: (value: number, offset: number) => dv.setInt32(offset, value, false),
            writeUInt16BE: (value: number, offset: number) => dv.setUint16(offset, value, false),
            writeUInt32BE: (value: number, offset: number) => dv.setUint32(offset, value, false),
            writeFloatBE: (value: number, offset: number) => dv.setFloat32(offset, value, false),
            writeDoubleBE: (value: number, offset: number) => dv.setFloat64(offset, value, false),

            copy: (targetBuffer: Uint8Array, targetStart?: number, sourceStart?: number, sourceEnd?: number) => {
                targetStart = defaults(targetStart, 0);
                sourceStart = defaults(sourceStart, 0);
                sourceEnd = defaults(sourceEnd, array.length);
                targetBuffer.set(array.subarray(sourceStart, sourceEnd), targetStart);
                return sourceEnd - sourceStart;
            }
        });
    }

    export function fromArrayBuffer(arrayBuffer: ArrayBuffer): SimpleBuffer {
        return fromUint8Array(new Uint8Array(arrayBuffer));
    }

    export function fromBuffer(buffer: Buffer): SimpleBuffer {
        return buffer;
    }

    export const IsNativeEndianLittle = new Uint16Array(new Uint8Array([0x12, 0x34]).buffer)[0] === 0x3412;

    /** source and target can't be the same */
    export function flipByteOrder(source: SimpleBuffer, target: Uint8Array, byteCount: number, elementByteSize: number, offset: number) {
        for (let i = 0, n = byteCount; i < n; i += elementByteSize) {
            for (let j = 0; j < elementByteSize; j++) {
                target[offset + i + elementByteSize - j - 1] = source[offset + i + j];
            }
        }
    }

    export function flipByteOrderInPlace2(buffer: ArrayBuffer, byteOffset = 0, length?: number) {
        const intView = new Int16Array(buffer, byteOffset, length);
        for (let i = 0, n = intView.length; i < n; ++i) {
            const val = intView[i];
            intView[i] = ((val & 0xff) << 8) | ((val >> 8) & 0xff);
        }
    }

    export function ensureLittleEndian(source: SimpleBuffer, target: SimpleBuffer, byteCount: number, elementByteSize: number, offset: number) {
        if (IsNativeEndianLittle) return;
        if (!byteCount || elementByteSize <= 1) return;
        flipByteOrder(source, target, byteCount, elementByteSize, offset);
    }
}