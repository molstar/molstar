/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from CIFTools.js (https://github.com/dsehnal/CIFTools.js)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../../../mol-data/db';
import StringBuilder from '../../../../mol-util/string-builder';
import { Category, Field, Encoder } from '../encoder';
import Writer from '../../writer';
import { getFieldDigitCount, getIncludedFields, getCategoryInstanceData, CategoryInstanceData } from './util';

export default class TextEncoder implements Encoder<string> {
    private builder = StringBuilder.create();
    private encoded = false;
    private dataBlockCreated = false;
    private filter: Category.Filter = Category.DefaultFilter;
    private formatter: Category.Formatter = Category.DefaultFormatter;

    readonly isBinary = false;

    binaryEncodingProvider = void 0;

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
        this.dataBlockCreated = true;
        StringBuilder.write(this.builder, `data_${(header || '').replace(/[ \n\t]/g, '').toUpperCase()}\n#\n`);
    }

    writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx, options?: Encoder.WriteCategoryOptions) {
        if (this.encoded) {
            throw new Error('The writer contents have already been encoded, no more writing.');
        }

        if (!this.dataBlockCreated) {
            throw new Error('No data block created.');
        }

        if (!options?.ignoreFilter && !this.filter.includeCategory(category.name)) return;
        const { instance, rowCount, source } = getCategoryInstanceData(category, context);
        if (!rowCount) return;

        if (rowCount === 1) {
            writeCifSingleRecord(category, instance, source, this.builder, this.filter, this.formatter);
        } else {
            writeCifLoop(category, instance, source, this.builder, this.filter, this.formatter);
        }
    }

    encode() {
        this.encoded = true;
    }

    writeTo(stream: Writer) {
        const chunks = StringBuilder.getChunks(this.builder);
        for (let i = 0, _i = chunks.length; i < _i; i++) {
            stream.writeString(chunks[i]);
        }
    }

    getSize() {
        return StringBuilder.getSize(this.builder);
    }

    getData() {
        return StringBuilder.getString(this.builder);
    }
}

function writeValue(builder: StringBuilder, data: any, key: any, f: Field<any, any>, floatPrecision: number, index: number): boolean {
    const kind = f.valueKind;
    const p = kind ? kind(key, data) : Column.ValueKind.Present;
    if (p !== Column.ValueKind.Present) {
        if (p === Column.ValueKind.NotPresent) writeNotPresent(builder);
        else writeUnknown(builder);
    } else {
        const val = f.value(key, data, index);
        const t = f.type;
        if (t === Field.Type.Str) {
            if (isMultiline(val as string)) {
                writeMultiline(builder, val as string);
                return true;
            } else {
                return writeChecked(builder, val as string);
            }
        } else if (t === Field.Type.Int) {
            writeInteger(builder, val as number);
        } else {
            writeFloat(builder, val as number, floatPrecision);
        }
    }
    return false;
}

function getFloatPrecisions(categoryName: string, fields: Field[], formatter: Category.Formatter) {
    const ret: number[] = [];
    for (const f of fields) {
        const format = formatter.getFormat(categoryName, f.name);
        if (format && typeof format.digitCount !== 'undefined') ret[ret.length] = f.type === Field.Type.Float ? Math.pow(10, Math.max(0, Math.min(format.digitCount, 15))) : 0;
        else ret[ret.length] = f.type === Field.Type.Float ? Math.pow(10, getFieldDigitCount(f)) : 0;
    }
    return ret;
}

function writeCifSingleRecord(category: Category, instance: Category.Instance, source: CategoryInstanceData['source'], builder: StringBuilder, filter: Category.Filter, formatter: Category.Formatter) {
    const fields = getIncludedFields(instance);
    const src = source[0];
    const data = src.data;
    let width = fields.reduce((w, f) => filter.includeField(category.name, f.name) ? Math.max(w, f.name.length) : 0, 0);

    // this means no field from this category is included.
    if (width === 0) return;
    width += category.name.length + 6;

    const it = src.keys();
    const key = it.move();
    const precisions = getFloatPrecisions(category.name, instance.fields, formatter);

    for (let _f = 0; _f < fields.length; _f++) {
        const f = fields[_f];
        if (!filter.includeField(category.name, f.name)) continue;

        StringBuilder.writePadRight(builder, `_${category.name}.${f.name}`, width);
        const multiline = writeValue(builder, data, key, f, precisions[_f], 0);
        if (!multiline) StringBuilder.newline(builder);
    }
    StringBuilder.write(builder, '#\n');
}

function writeCifLoop(category: Category, instance: Category.Instance, source: CategoryInstanceData['source'], builder: StringBuilder, filter: Category.Filter, formatter: Category.Formatter) {
    const fieldSource = getIncludedFields(instance);
    const fields = filter === Category.DefaultFilter ? fieldSource : fieldSource.filter(f => filter.includeField(category.name, f.name));
    const fieldCount = fields.length;
    if (fieldCount === 0) return;

    const precisions = getFloatPrecisions(category.name, fields, formatter);

    writeLine(builder, 'loop_');
    for (let i = 0; i < fieldCount; i++) {
        writeLine(builder, `_${category.name}.${fields[i].name}`);
    }

    let index = 0;
    for (let _c = 0; _c < source.length; _c++) {
        const src = source[_c];
        const data = src.data;

        if (src.rowCount === 0) continue;

        const it = src.keys();
        while (it.hasNext)  {
            const key = it.move();

            let multiline = false;
            for (let _f = 0; _f < fieldCount; _f++) {
                multiline = writeValue(builder, data, key, fields[_f], precisions[_f], index);
            }
            if (!multiline) StringBuilder.newline(builder);
            index++;
        }
    }
    StringBuilder.write(builder, '#\n');
}

function isMultiline(value: string) {
    return typeof value === 'string' && value.indexOf('\n') >= 0;
}

function writeLine(builder: StringBuilder, val: string) {
    StringBuilder.write(builder, val);
    StringBuilder.newline(builder);
}

function writeInteger(builder: StringBuilder, val: number) {
    StringBuilder.writeInteger(builder, val);
    StringBuilder.whitespace1(builder);
}

function writeFloat(builder: StringBuilder, val: number, precisionMultiplier: number) {
    StringBuilder.writeFloat(builder, val, precisionMultiplier);
    StringBuilder.whitespace1(builder);
}

function writeNotPresent(builder: StringBuilder) {
    StringBuilder.writeSafe(builder, '. ');
}

function writeUnknown(builder: StringBuilder) {
    StringBuilder.writeSafe(builder, '? ');
}

function writeChecked(builder: StringBuilder, val: string) {
    if (!val) {
        StringBuilder.writeSafe(builder, '. ');
        return false;
    }

    let escape = val.charCodeAt(0) === 95 /* _ */, escapeCharStart = '\'', escapeCharEnd = '\' ';
    let hasWhitespace = false;
    let hasSingle = false;
    let hasDouble = false;
    for (let i = 0, _l = val.length - 1; i < _l; i++) {
        const c = val.charCodeAt(i);

        switch (c) {
            case 9: hasWhitespace = true; break; // \t
            case 10: // \n
                writeMultiline(builder, val);
                return true;
            case 32: hasWhitespace = true; break; // ' '
            case 34: // "
                if (hasSingle) {
                    writeMultiline(builder, val);
                    return true;
                }

                hasDouble = true;
                escape = true;
                escapeCharStart = '\'';
                escapeCharEnd = '\' ';
                break;
            case 39: // '
                if (hasDouble) {
                    writeMultiline(builder, val);
                    return true;
                }

                escape = true;
                hasSingle = true;
                escapeCharStart = '"';
                escapeCharEnd = '" ';
                break;
        }
    }

    const fst = val.charCodeAt(0);
    if (!escape && (fst === 35 /* # */|| fst === 36 /* $ */ || fst === 59 /* ; */ || fst === 91 /* [ */ || fst === 93 /* ] */ || hasWhitespace)) {
        escapeCharStart = '\'';
        escapeCharEnd = '\' ';
        escape = true;
    }

    if (escape) {
        StringBuilder.writeSafe(builder, escapeCharStart);
        StringBuilder.writeSafe(builder, val);
        StringBuilder.writeSafe(builder, escapeCharEnd);
    } else {
        StringBuilder.writeSafe(builder, val);
        StringBuilder.writeSafe(builder, ' ');
    }

    return false;
}

function writeMultiline(builder: StringBuilder, val: string) {
    StringBuilder.writeSafe(builder, '\n;' + val);
    StringBuilder.writeSafe(builder, '\n;\n');
}
