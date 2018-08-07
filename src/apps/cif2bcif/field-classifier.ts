/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-data/db'
import { CifField } from 'mol-io/reader/cif/data-model'
import { CifWriter } from 'mol-io/writer/cif'
import { ArrayEncoder, ArrayEncoding as E } from 'mol-io/common/binary-cif';

namespace IntClassifier {
    function packSize(value: number, upperLimit: number) {
        return value >= 0
            ? Math.ceil((value + 1) / upperLimit)
            : Math.ceil((value + 1) / (-upperLimit - 1));
    }

    type IntColumnInfo = { signed: boolean, limit8: number, limit16: number };

    function getInfo(data: number[]): IntColumnInfo {
        let signed = false;
        for (let i = 0, n = data.length; i < n; i++) {
            if (data[i] < 0) {
                signed = true;
                break;
            }
        }
        return signed ? { signed, limit8: 0x7F, limit16: 0x7FFF } : { signed, limit8: 0xFF, limit16: 0xFFFF };
    }

    type SizeInfo = { pack8: number, pack16: number, count: number }
    function SizeInfo(): SizeInfo { return { pack8: 0, pack16: 0, count: 0 } };

    function incSize({ limit8, limit16 }: IntColumnInfo, info: SizeInfo, value: number) {
        info.pack8 += packSize(value, limit8);
        info.pack16 += packSize(value, limit16);
        info.count += 1;
    }

    function incSizeSigned(info: SizeInfo, value: number) {
        info.pack8 += packSize(value, 0x7F);
        info.pack16 += packSize(value, 0x7FFF);
        info.count += 1;
    }

    function byteSize(info: SizeInfo) {
        if (info.count * 4 < info.pack16 * 2) return { length: info.count * 4, elem: 4 };
        if (info.pack16 * 2 < info.pack8) return { length: info.pack16 * 2, elem: 2 };
        return { length: info.pack8, elem: 1 };
    }

    function packingSize(data: number[], info: IntColumnInfo) {
        const size = SizeInfo();
        for (let i = 0, n = data.length; i < n; i++) {
            incSize(info, size, data[i]);
        }
        return { ...byteSize(size), kind: 'pack' };
    }

    function deltaSize(data: number[], info: IntColumnInfo) {
        const size = SizeInfo();
        let prev = data[0];
        for (let i = 1, n = data.length; i < n; i++) {
            incSizeSigned(size, data[i] - prev);
            prev = data[i];
        }
        return { ...byteSize(size), kind: 'delta' };
    }

    function rleSize(data: number[], info: IntColumnInfo) {
        const size = SizeInfo();
        let run = 1;
        for (let i = 1, n = data.length; i < n; i++) {
            if (data[i - 1] !== data[i]) {
                incSize(info, size, data[i - 1]);
                incSize(info, size, run);
                run = 1;
            } else {
                run++;
            }
        }
        incSize(info, size, data[data.length - 1]);
        incSize(info, size, run);

        return { ...byteSize(size), kind: 'rle' };
    }

    function deltaRleSize(data: number[], info: IntColumnInfo) {
        const size = SizeInfo();
        let run = 1, prev = 0, prevValue = 0;
        for (let i = 1, n = data.length; i < n; i++) {
            const v = data[i] - prev;
            if (prevValue !== v) {
                incSizeSigned(size, prevValue);
                incSizeSigned(size, run);
                run = 1;
            } else {
                run++;
            }
            prevValue = v;
            prev = data[i];
        }
        incSizeSigned(size, prevValue);
        incSizeSigned(size, run);

        return { ...byteSize(size), kind: 'delta-rle' };
    }

    export function getSize(data: number[]) {
        const info = getInfo(data);
        const sizes = [packingSize(data, info), rleSize(data, info), deltaSize(data, info), deltaRleSize(data, info)];
        sizes.sort((a, b) => a.length - b.length);
        return sizes;
    }

    export function classify(data: number[], name: string): ArrayEncoder {
        if (data.length < 2) return E.by(E.byteArray);

        const sizes = getSize(data);
        const size = sizes[0];
        // console.log(`${name}: ${size.kind} ${size.length}b ${data.length}`);
        // console.log(`${name}: ${sizes.map(s => `${s.kind}: ${s.length}b`).join(' | ')}`);

        switch (size.kind) {
            case 'pack': return E.by(E.integerPacking);
            case 'rle': return E.by(E.runLength).and(E.integerPacking);
            case 'delta': return E.by(E.delta).and(E.integerPacking);
            case 'delta-rle': return E.by(E.delta).and(E.runLength).and(E.integerPacking);
        }

        throw 'bug';
    }
}

namespace FloatClassifier {
    const delta = 1e-6;
    function digitCount(v: number) {
        let m = 1;
        for (let i = 0; i < 5; i++) {
            const r = Math.round(m * v) / m;
            if (Math.abs(v - r) < delta) return m;
            m *= 10;
        }
        return 10000;
    }

    export function classify(data: number[], name: string) {
        // if a vector/matrix, do not reduce precision
        if (name.indexOf('[') > 0) return { encoder: E.by(E.byteArray), typedArray: Float64Array };

        let dc = 10;
        for (let i = 0, n = data.length; i < n; i++) dc = Math.max(dc, digitCount(data[i]));

        if (dc >= 10000) return { encoder: E.by(E.byteArray), typedArray: Float64Array };

        const intArray = new Int32Array(data.length);
        for (let i = 0, n = data.length; i < n; i++) intArray[i] = data[i] * dc;

        const sizes = IntClassifier.getSize(intArray as any);
        const size = sizes[0];

        // console.log(`>> ${name}: ${size.kind} ${size.length}b ${data.length} x${dc}`);
        // console.log(`   ${name}: ${sizes.map(s => `${s.kind}: ${s.length}b`).join(' | ')}`);

        switch (size.kind) {
            case 'pack': return { encoder: E.by(E.fixedPoint(dc)).and(E.integerPacking), typedArray: Float32Array };
            case 'rle': return { encoder: E.by(E.fixedPoint(dc)).and(E.runLength).and(E.integerPacking), typedArray: Float32Array };
            case 'delta': return { encoder: E.by(E.fixedPoint(dc)).and(E.delta).and(E.integerPacking), typedArray: Float32Array };
            case 'delta-rle': return { encoder: E.by(E.fixedPoint(dc)).and(E.delta).and(E.runLength).and(E.integerPacking), typedArray: Float32Array };
        }

        throw 'bug';
    }
}

const intRegex = /^-?\d+$/
const floatRegex = /^-?(([0-9]+)[.]?|([0-9]*[.][0-9]+))([(][0-9]+[)])?([eE][+-]?[0-9]+)?$/

// Classify a cif field as str, int or float based the data it contains.
// To classify a field as int or float all items are checked.
function classify(name: string, field: CifField): CifWriter.Field {
    let floatCount = 0, hasString = false;
    for (let i = 0, _i = field.rowCount; i < _i; i++) {
        const k = field.valueKind(i);
        if (k !== Column.ValueKind.Present) continue;
        const v = field.str(i);
        if (intRegex.test(v)) continue;
        else if (floatRegex.test(v)) floatCount++;
        else { hasString = true; break; }
    }

    if (hasString) return { name, type: CifWriter.Field.Type.Str, value: field.str, valueKind: field.valueKind };
    if (floatCount > 0) {
        const { encoder, typedArray } = FloatClassifier.classify(field.toFloatArray({ array: Float64Array }) as number[], name)
        return CifWriter.Field.float(name, field.float, { valueKind: field.valueKind, encoder, typedArray });
    } else {
        const encoder = IntClassifier.classify(field.toIntArray({ array: Int32Array }) as number[], name);
        return CifWriter.Field.int(name, field.int, { valueKind: field.valueKind, encoder, typedArray: Int32Array });
    }
}

export default classify;