/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from CIFTools.js (https://github.com/dsehnal/CIFTools.js; MIT) and MMTF (https://github.com/rcsb/mmtf-javascript/; MIT)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import ChunkedArray from 'mol-base/collections/chunked-array'
import { Encoding, EncodedData } from './encoding'

export interface ArrayEncoder {
    and(f: ArrayEncoding.Provider): ArrayEncoder,
    encode(data: ArrayLike<any>): EncodedData
}

export class ArrayEncoderImpl implements ArrayEncoder {
    and(f: ArrayEncoding.Provider) {
        return new ArrayEncoderImpl(this.providers.concat([f]));
    }

    encode(data: ArrayLike<any>): EncodedData {
        let encoding: Encoding[] = [];
        for (let p of this.providers) {
            let t = p(data);

            if (!t.encodings.length) {
                throw new Error('Encodings must be non-empty.');
            }

            data = t.data;
            for (let e of t.encodings) {
                encoding.push(e);
            }
        }
        if (!(data instanceof Uint8Array)) {
            throw new Error('The encoding must result in a Uint8Array. Fix your encoding chain.');
        }
        return {
            encoding,
            data
        }
    }

    constructor(private providers: ArrayEncoding.Provider[]) {

    }
}

export namespace ArrayEncoder {
    export function by(f: ArrayEncoding.Provider): ArrayEncoder {
        return new ArrayEncoderImpl([f]);
    }
}

export namespace ArrayEncoding {
    export type TypedArrayCtor = { new(size: number): ArrayLike<number> & { buffer: ArrayBuffer, byteLength: number, byteOffset: number, BYTES_PER_ELEMENT: number } }

    export interface Result {
        encodings: Encoding[],
        data: any
    }

    export type Provider = (data: any) => Result

    export function by(f: Provider): ArrayEncoder {
        return new ArrayEncoderImpl([f]);
    }

    function uint8(data: Uint8Array): Result {
        return {
            encodings: [{ kind: 'ByteArray', type: Encoding.IntDataType.Uint8 }],
            data
        };
    }

    function int8(data: Int8Array): Result {
        return {
            encodings: [{ kind: 'ByteArray', type: Encoding.IntDataType.Int8 }],
            data: new Uint8Array(data.buffer, data.byteOffset)
        };
    }

    const writers = {
        [Encoding.IntDataType.Int16]: function (v: DataView, i: number, a: number) { v.setInt16(2 * i, a, true) },
        [Encoding.IntDataType.Uint16]: function (v: DataView, i: number, a: number) { v.setUint16(2 * i, a, true) },
        [Encoding.IntDataType.Int32]: function (v: DataView, i: number, a: number) { v.setInt32(4 * i, a, true) },
        [Encoding.IntDataType.Uint32]: function (v: DataView, i: number, a: number) { v.setUint32(4 * i, a, true) },
        [Encoding.FloatDataType.Float32]: function (v: DataView, i: number, a: number) { v.setFloat32(4 * i, a, true) },
        [Encoding.FloatDataType.Float64]: function (v: DataView, i: number, a: number) { v.setFloat64(8 * i, a, true) }
    }

    const byteSizes = {
        [Encoding.IntDataType.Int16]: 2,
        [Encoding.IntDataType.Uint16]: 2,
        [Encoding.IntDataType.Int32]: 4,
        [Encoding.IntDataType.Uint32]: 4,
        [Encoding.FloatDataType.Float32]: 4,
        [Encoding.FloatDataType.Float64]: 8
    }

    export function byteArray(data: Encoding.FloatArray | Encoding.IntArray) {
        let type = Encoding.getDataType(data);

        if (type === Encoding.IntDataType.Int8) return int8(data as Int8Array);
        else if (type === Encoding.IntDataType.Uint8) return uint8(data as Uint8Array);

        let result = new Uint8Array(data.length * byteSizes[type]);
        let w = writers[type];
        let view = new DataView(result.buffer);
        for (let i = 0, n = data.length; i < n; i++) {
            w(view, i, data[i]);
        }
        return {
            encodings: [<Encoding.ByteArray>{ kind: 'ByteArray', type }],
            data: result
        };
    }

    function _fixedPoint(data: Encoding.FloatArray, factor: number): Result {
        let srcType = Encoding.getDataType(data) as Encoding.FloatDataType;
        let result = new Int32Array(data.length);
        for (let i = 0, n = data.length; i < n; i++) {
            result[i] = Math.round(data[i] * factor);
        }
        return {
            encodings: [{ kind: 'FixedPoint', factor, srcType }],
            data: result
        };
    }
    export function fixedPoint(factor: number): Provider { return data => _fixedPoint(data as Encoding.FloatArray, factor); }

    function _intervalQuantizaiton(data: Encoding.FloatArray, min: number, max: number, numSteps: number, arrayType: new (size: number) => Encoding.IntArray): Result {
        let srcType = Encoding.getDataType(data) as Encoding.FloatDataType;
        if (!data.length) {
            return {
                encodings: [{ kind: 'IntervalQuantization', min, max, numSteps, srcType }],
                data: new Int32Array(0)
            };
        }

        if (max < min) {
            let t = min;
            min = max;
            max = t;
        }

        let delta = (max - min) / (numSteps - 1);

        let output = new arrayType(data.length);
        for (let i = 0, n = data.length; i < n; i++) {
            let v = data[i];
            if (v <= min) output[i] = 0;
            else if (v >= max) output[i] = numSteps;
            else output[i] = (Math.round((v - min) / delta)) | 0;
        }

        return {
            encodings: [{ kind: 'IntervalQuantization', min, max, numSteps, srcType }],
            data: output
        };
    }
    export function intervalQuantizaiton(min: number, max: number, numSteps: number, arrayType: new (size: number) => Encoding.IntArray = Int32Array): Provider {
        return data => _intervalQuantizaiton(data as Encoding.FloatArray, min, max, numSteps, arrayType);
    }

    export function runLength(data: Encoding.IntArray): Result {
        let srcType = Encoding.getDataType(data) as Encoding.IntDataType;
        if (srcType === void 0) {
            data = new Int32Array(data);
            srcType = Encoding.IntDataType.Int32;
        }

        if (!data.length) {
            return {
                encodings: [{ kind: 'RunLength', srcType, srcSize: 0 }],
                data: new Int32Array(0)
            };
        }

        // calculate output size
        let fullLength = 2;
        for (let i = 1, il = data.length; i < il; i++) {
            if (data[i - 1] !== data[i]) {
                fullLength += 2;
            }
        }
        let output = new Int32Array(fullLength);
        let offset = 0;
        let runLength = 1;
        for (let i = 1, il = data.length; i < il; i++) {
            if (data[i - 1] !== data[i]) {
                output[offset] = data[i - 1];
                output[offset + 1] = runLength;
                runLength = 1;
                offset += 2;
            } else {
                ++runLength;
            }
        }
        output[offset] = data[data.length - 1];
        output[offset + 1] = runLength;
        return {
            encodings: [{ kind: 'RunLength', srcType, srcSize: data.length }],
            data: output
        };
    }

    export function delta(data: Int8Array | Int16Array | Int32Array): Result {
        if (!Encoding.isSignedIntegerDataType(data)) {
            throw new Error('Only signed integer types can be encoded using delta encoding.');
        }

        let srcType = Encoding.getDataType(data) as Encoding.IntDataType;
        if (srcType === void 0) {
            data = new Int32Array(data);
            srcType = Encoding.IntDataType.Int32;
        }
        if (!data.length) {
            return {
                encodings: [{ kind: 'Delta', origin: 0, srcType }],
                data: new (data as any).constructor(0)
            };
        }

        let output = new (data as any).constructor(data.length);
        let origin = data[0];
        output[0] = data[0];
        for (let i = 1, n = data.length; i < n; i++) {
            output[i] = data[i] - data[i - 1];
        }
        output[0] = 0;
        return {
            encodings: [{ kind: 'Delta', origin, srcType }],
            data: output
        };
    }

    function isSigned(data: Int32Array) {
        for (let i = 0, n = data.length; i < n; i++) {
            if (data[i] < 0) return true;
        }
        return false;
    }

    function packingSize(data: Int32Array, upperLimit: number) {
        let lowerLimit = -upperLimit - 1;
        let size = 0;
        for (let i = 0, n = data.length; i < n; i++) {
            let value = data[i];
            if (value === 0) {
                size += 1;
            } else if (value > 0) {
                size += Math.ceil(value / upperLimit);
                if (value % upperLimit === 0) size += 1;
            } else {
                size += Math.ceil(value / lowerLimit);
                if (value % lowerLimit === 0) size += 1;
            }
        }
        return size;
    }

    function determinePacking(data: Int32Array): { isSigned: boolean, size: number, bytesPerElement: number } {
        let signed = isSigned(data);
        let size8 = signed ? packingSize(data, 0x7F) : packingSize(data, 0xFF);
        let size16 = signed ? packingSize(data, 0x7FFF) : packingSize(data, 0xFFFF);

        if (data.length * 4 < size16 * 2) {
            // 4 byte packing is the most effective
            return {
                isSigned: signed,
                size: data.length,
                bytesPerElement: 4
            };
        } else if (size16 * 2 < size8) {
            // 2 byte packing is the most effective
            return {
                isSigned: signed,
                size: size16,
                bytesPerElement: 2
            }
        } else {
            // 1 byte packing is the most effective
            return {
                isSigned: signed,
                size: size8,
                bytesPerElement: 1
            }
        };
    }

    function _integerPacking(data: Int32Array, packing: { isSigned: boolean, size: number, bytesPerElement: number }): Result {
        let upperLimit = packing.isSigned
            ? (packing.bytesPerElement === 1 ? 0x7F : 0x7FFF)
            : (packing.bytesPerElement === 1 ? 0xFF : 0xFFFF);

        let lowerLimit = -upperLimit - 1;
        let n = data.length;
        let packed = packing.isSigned
            ? packing.bytesPerElement === 1 ? new Int8Array(packing.size) : new Int16Array(packing.size)
            : packing.bytesPerElement === 1 ? new Uint8Array(packing.size) : new Uint16Array(packing.size);
        let j = 0;
        for (let i = 0; i < n; i++) {
            let value = data[i];
            if (value >= 0) {
                while (value >= upperLimit) {
                    packed[j] = upperLimit;
                    ++j;
                    value -= upperLimit;
                }
            } else {
                while (value <= lowerLimit) {
                    packed[j] = lowerLimit;
                    ++j;
                    value -= lowerLimit;
                }
            }
            packed[j] = value;
            ++j;
        }

        let result = byteArray(packed);
        return {
            encodings: [{
                kind: 'IntegerPacking',
                byteCount: packing.bytesPerElement,
                isUnsigned: !packing.isSigned,
                srcSize: n
            },
            result.encodings[0]
            ],
            data: result.data
        };
    }

    /**
     * Packs Int32 array. The packing level is determined automatically to either 1-, 2-, or 4-byte words.
     */
    export function integerPacking(data: Int32Array): Result {
        if (!(data instanceof Int32Array)) {
            throw new Error('Integer packing can only be applied to Int32 data.');
        }

        let packing = determinePacking(data);

        if (packing.bytesPerElement === 4) {
            // no packing done, Int32 encoding will be used
            return byteArray(data);
        }

        return _integerPacking(data, packing);
    }

    export function stringArray(data: string[]): Result {
        let map: any = Object.create(null);
        let strings: string[] = [];
        let accLength = 0;
        let offsets = ChunkedArray.create<number>(s => new Int32Array(s), 1, 1024, true)
        let output = new Int32Array(data.length);

        ChunkedArray.add(offsets, 0);
        let i = 0;
        for (let s of data) {
            // handle null strings.
            if (s === null || s === void 0) {
                output[i++] = -1;
                continue;
            }

            let index = map[s];
            if (index === void 0) {
                // increment the length
                accLength += s.length;

                // store the string and index
                index = strings.length;
                strings[index] = s;
                map[s] = index;

                // write the offset
                ChunkedArray.add(offsets, accLength);
            }
            output[i++] = index;
        }

        let encOffsets = ArrayEncoder.by(delta).and(integerPacking).encode(ChunkedArray.compact(offsets));
        let encOutput = ArrayEncoder.by(delta).and(runLength).and(integerPacking).encode(output);

        return {
            encodings: [{ kind: 'StringArray', dataEncoding: encOutput.encoding, stringData: strings.join(''), offsetEncoding: encOffsets.encoding, offsets: encOffsets.data }],
            data: encOutput.data
        };
    }
}