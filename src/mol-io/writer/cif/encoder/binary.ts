/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from CIFTools.js (https://github.com/dsehnal/CIFTools.js)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../../../mol-data/db';
import encodeMsgPack from '../../../common/msgpack/encode';
import {
    EncodedColumn, EncodedData, EncodedFile, EncodedDataBlock, EncodedCategory, ArrayEncoder, ArrayEncoding as E, VERSION
} from '../../../common/binary-cif';
import { Field, Category, Encoder } from '../encoder';
import Writer from '../../writer';
import { getIncludedFields, getCategoryInstanceData, CategoryInstanceData } from './util';
import { classifyIntArray, classifyFloatArray } from '../../../common/binary-cif/classifier';
import { ArrayCtor } from '../../../../mol-util/type-helpers';

export interface BinaryEncodingProvider {
    get(category: string, field: string): ArrayEncoder | undefined;
}

export default class BinaryEncoder implements Encoder<Uint8Array> {
    private data: EncodedFile;
    private dataBlocks: EncodedDataBlock[] = [];
    private encodedData: Uint8Array;
    private filter: Category.Filter = Category.DefaultFilter;
    private formatter: Category.Formatter = Category.DefaultFormatter;

    readonly isBinary = true;

    binaryEncodingProvider: BinaryEncodingProvider | undefined = void 0;

    setFilter(filter?: Category.Filter) {
        this.filter = filter || Category.DefaultFilter;
    }

    isCategoryIncluded(name: string) {
        return this.filter.includeCategory(name);
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

    writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx, options?: Encoder.WriteCategoryOptions) {
        if (!this.data) {
            throw new Error('The writer contents have already been encoded, no more writing.');
        }

        if (!this.dataBlocks.length) {
            throw new Error('No data block created.');
        }

        if (!options?.ignoreFilter && !this.filter.includeCategory(category.name)) return;

        const { instance, rowCount, source } = getCategoryInstanceData(category, context);
        if (!rowCount) return;

        const cat: EncodedCategory = { name: '_' + category.name, columns: [], rowCount };
        const fields = getIncludedFields(instance);

        for (const f of fields) {
            if (!this.filter.includeField(category.name, f.name)) continue;

            const format = this.formatter.getFormat(category.name, f.name);
            cat.columns.push(encodeField(category.name, f, source, rowCount, format, this.binaryEncodingProvider, this.autoClassify));
        }
        // no columns included.
        if (!cat.columns.length) return;

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

    getSize() {
        return this.encodedData.length;
    }

    constructor(encoder: string, encodingProvider: BinaryEncodingProvider | undefined, private autoClassify: boolean) {
        this.binaryEncodingProvider = encodingProvider;
        this.data = {
            encoder,
            version: VERSION,
            dataBlocks: this.dataBlocks
        };
    }
}

function getArrayCtor(field: Field, format: Field.Format | undefined): ArrayCtor<string | number> {
    if (format && format.typedArray) return format.typedArray;
    if (field.defaultFormat && field.defaultFormat.typedArray) return field.defaultFormat.typedArray;
    if (field.type === Field.Type.Str) return Array;
    if (field.type === Field.Type.Int) return Int32Array;
    return Float64Array;
}

function getDefaultEncoder(type: Field.Type): ArrayEncoder {
    if (type === Field.Type.Str) return ArrayEncoder.by(E.stringArray);
    return ArrayEncoder.by(E.byteArray);
}

function tryGetEncoder(categoryName: string, field: Field, format: Field.Format | undefined, provider: BinaryEncodingProvider | undefined) {
    if (format && format.encoder) {
        return format.encoder;
    } else if (field.defaultFormat && field.defaultFormat.encoder) {
        return field.defaultFormat.encoder;
    } else if (provider) {
        return provider.get(categoryName, field.name);
    } else {
        return void 0;
    }
}

function classify(type: Field.Type, data: ArrayLike<any>) {
    if (type === Field.Type.Str) return ArrayEncoder.by(E.stringArray);
    if (type === Field.Type.Int) return classifyIntArray(data);
    return classifyFloatArray(data);
}

function encodeField(categoryName: string, field: Field, data: CategoryInstanceData['source'], totalCount: number,
    format: Field.Format | undefined, encoderProvider: BinaryEncodingProvider | undefined, autoClassify: boolean): EncodedColumn {

    const { array, allPresent, mask } = getFieldData(field, getArrayCtor(field, format), totalCount, data);

    let encoder: ArrayEncoder | undefined = tryGetEncoder(categoryName, field, format, encoderProvider);
    if (!encoder) {
        if (autoClassify) encoder = classify(field.type, array);
        else encoder = getDefaultEncoder(field.type);
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

function getFieldData(field: Field<any, any>, arrayCtor: ArrayCtor<string | number>, totalCount: number, data: CategoryInstanceData['source']) {
    const isStr = field.type === Field.Type.Str;
    const array = new arrayCtor(totalCount);
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
                if (isStr)
                    array[offset] = '';
                allPresent = false;
            } else {
                mask[offset] = Column.ValueKind.Present;
                array[offset] = getter(key, d, offset);
            }
            offset++;
        }
    }
    return { array, allPresent, mask };
}
