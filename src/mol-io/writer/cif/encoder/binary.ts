/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from CIFTools.js (https://github.com/dsehnal/CIFTools.js)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Iterator } from 'mol-data'
import { Column } from 'mol-data/db'
import encodeMsgPack from '../../../common/msgpack/encode'
import {
    EncodedColumn, EncodedData, EncodedFile, EncodedDataBlock, EncodedCategory, ArrayEncoder, ArrayEncoding as E, VERSION, ArrayEncoding
} from '../../../common/binary-cif'
import { Field, Category, Encoder } from '../encoder'
import Writer from '../../writer'

export default class BinaryEncoder implements Encoder<Uint8Array> {
    private data: EncodedFile;
    private dataBlocks: EncodedDataBlock[] = [];
    private encodedData: Uint8Array;
    private filter: Category.Filter = Category.DefaultFilter;
    private formatter: Category.Formatter = Category.DefaultFormatter;

    setFilter(filter?: Category.Filter) {
        this.filter = filter || Category.DefaultFilter;
    }

    setFormatter(formatter?: Category.Formatter) {
        this.formatter = formatter || Category.DefaultFormatter;
    }

    startDataBlock(header: string) {
        this.dataBlocks.push({
            header: (header || '').replace(/[ \n\t]/g, '').toUpperCase(),
            categories: []
        });
    }

    writeCategory<Ctx>(category: Category.Provider<Ctx>, contexts?: Ctx[]) {
        if (!this.data) {
            throw new Error('The writer contents have already been encoded, no more writing.');
        }

        if (!this.dataBlocks.length) {
            throw new Error('No data block created.');
        }

        const src = !contexts || !contexts.length ? [category(<any>void 0)] : contexts.map(c => category(c));
        const categories = src.filter(c => c && c.rowCount > 0);
        if (!categories.length) return;
        if (!this.filter.includeCategory(categories[0].name)) return;

        const count = categories.reduce((a, c) => a + c.rowCount, 0);
        if (!count) return;

        const first = categories[0]!;
        const cat: EncodedCategory = { name: '_' + first.name, columns: [], rowCount: count };
        const data = categories.map(c => ({ data: c.data, keys: () => c.keys ? c.keys() : Iterator.Range(0, c.rowCount - 1) }));
        for (const f of first.fields) {
            if (!this.filter.includeField(first.name, f.name)) continue;

            const format = this.formatter.getFormat(first.name, f.name);
            cat.columns.push(encodeField(f, data, count, getArrayCtor(f, format), getEncoder(f, format)));
        }
        this.dataBlocks[this.dataBlocks.length - 1].categories.push(cat);
    }

    encode() {
        if (this.encodedData) return;
        this.encodedData = encodeMsgPack(this.data);
        this.data = <any>null;
        this.dataBlocks = <any>null;
    }

    writeTo(writer: Writer) {
        writer.writeBinary(this.encodedData);
    }

    getData() {
        this.encode();
        return this.encodedData;
    }

    constructor(encoder: string) {
        this.data = {
            encoder,
            version: VERSION,
            dataBlocks: this.dataBlocks
        };
    }
}

function createArray(type: Field.Type, arrayCtor: ArrayEncoding.TypedArrayCtor | undefined,  count: number) {
    if (type === Field.Type.Str) return new Array(count) as any;
    else if (arrayCtor) return new arrayCtor(count) as any;
    else return (type === Field.Type.Int ? new Int32Array(count) : new Float32Array(count)) as any;
}

function getArrayCtor(field: Field, format: Field.Format | undefined) {
    if (format && format.typedArray) return format.typedArray;
    if (field.defaultFormat && field.defaultFormat.typedArray) return field.defaultFormat.typedArray;
    return void 0;
}

function getEncoder(field: Field, format: Field.Format | undefined) {
    if (format && format.encoder) return format.encoder;
    if (field.defaultFormat && field.defaultFormat.encoder) {
        return field.defaultFormat.encoder;
    } else if (field.type === Field.Type.Str) {
        return ArrayEncoder.by(E.stringArray);
    } else {
        return ArrayEncoder.by(E.byteArray);
    }
}

function encodeField(field: Field, data: { data: any, keys: () => Iterator<any> }[], totalCount: number, arrayCtor: ArrayEncoding.TypedArrayCtor | undefined, encoder: ArrayEncoder): EncodedColumn {
    const isStr = field.type === Field.Type.Str;
    const array = createArray(field.type, arrayCtor, totalCount);
    const mask = new Uint8Array(totalCount);
    const valueKind = field.valueKind;
    const getter = field.value;
    let allPresent = true;

    let offset = 0;
    for (let _d = 0; _d < data.length; _d++) {
        const d = data[_d].data;
        const keys = data[_d].keys();
        while (keys.hasNext) {
            const key = keys.move();
            const p = valueKind ? valueKind(key, d) : Column.ValueKind.Present;
            if (p !== Column.ValueKind.Present) {
                mask[offset] = p;
                if (isStr) array[offset] = '';
                allPresent = false;
            } else {
                mask[offset] = Column.ValueKind.Present;
                array[offset] = getter(key, d);
            }
            offset++;
        }
    }

    const encoded = encoder.encode(array);

    let maskData: EncodedData | undefined = void 0;

    if (!allPresent) {
        const maskRLE = ArrayEncoder.by(E.runLength).and(E.byteArray).encode(mask);
        if (maskRLE.data.length < mask.length) {
            maskData = maskRLE;
        } else {
            maskData = ArrayEncoder.by(E.byteArray).encode(mask);
        }
    }

    return {
        name: field.name,
        data: encoded,
        mask: maskData
    };
}