/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from CIFTools.js (https://github.com/dsehnal/CIFTools.js)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Iterator } from 'mol-data'
import { Column } from 'mol-data/db'
import StringBuilder from 'mol-util/string-builder'
import { Category, Field, Encoder } from '../encoder'
import Writer from '../../writer'
import { getFieldDigitCount, getIncludedFields } from './util';

export default class TextEncoder implements Encoder<string> {
    private builder = StringBuilder.create();
    private encoded = false;
    private dataBlockCreated = false;
    private filter: Category.Filter = Category.DefaultFilter;
    private formatter: Category.Formatter = Category.DefaultFormatter;

    setFilter(filter?: Category.Filter) {
        this.filter = filter || Category.DefaultFilter;
    }

    setFormatter(formatter?: Category.Formatter) {
        this.formatter = formatter || Category.DefaultFormatter;
    }

    startDataBlock(header: string) {
        this.dataBlockCreated = true;
        StringBuilder.write(this.builder, `data_${(header || '').replace(/[ \n\t]/g, '').toUpperCase()}\n#\n`);
    }

    writeCategory<Ctx>(category: Category.Provider<Ctx>, contexts?: Ctx[]) {
        if (this.encoded) {
            throw new Error('The writer contents have already been encoded, no more writing.');
        }

        if (!this.dataBlockCreated) {
            throw new Error('No data block created.');
        }

        const categories = !contexts || !contexts.length ? [category(<any>void 0)] : contexts.map(c => category(c));
        if (!categories.length) return;
        if (!this.filter.includeCategory(categories[0].name)) return;

        const rowCount = categories.reduce((v, c) => v + c.rowCount, 0);

        if (rowCount === 0) return;

        if (rowCount === 1) {
            writeCifSingleRecord(categories[0]!, this.builder, this.filter, this.formatter);
        } else {
            writeCifLoop(categories, this.builder, this.filter, this.formatter);
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

    getData() {
        return StringBuilder.getString(this.builder);
    }
}

function writeValue(builder: StringBuilder, data: any, key: any, f: Field<any, any>, floatPrecision: number): boolean {
    const kind = f.valueKind;
    const p = kind ? kind(key, data) : Column.ValueKind.Present;
    if (p !== Column.ValueKind.Present) {
        if (p === Column.ValueKind.NotPresent) writeNotPresent(builder);
        else writeUnknown(builder);
    } else {
        const val = f.value(key, data);
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

function writeCifSingleRecord(category: Category<any>, builder: StringBuilder, filter: Category.Filter, formatter: Category.Formatter) {
    const fields = getIncludedFields(category);
    const data = category.data;
    let width = fields.reduce((w, f) => filter.includeField(category.name, f.name) ? Math.max(w, f.name.length) : 0, 0);

    // this means no field from this category is included.
    if (width === 0) return;
    width += category.name.length + 6;

    const it = category.keys ? category.keys() : Iterator.Range(0, category.rowCount - 1);
    const key = it.move();

    const precisions = getFloatPrecisions(category.name, category.fields, formatter);

    for (let _f = 0; _f < fields.length; _f++) {
        const f = fields[_f];
        if (!filter.includeField(category.name, f.name)) continue;

        StringBuilder.writePadRight(builder, `_${category.name}.${f.name}`, width);
        const multiline = writeValue(builder, data, key, f, precisions[_f]);
        if (!multiline) StringBuilder.newline(builder);
    }
    StringBuilder.write(builder, '#\n');
}

function writeCifLoop(categories: Category[], builder: StringBuilder, filter: Category.Filter, formatter: Category.Formatter) {
    const first = categories[0];
    const fieldSource = getIncludedFields(first);
    const fields = filter === Category.DefaultFilter ? fieldSource : fieldSource.filter(f => filter.includeField(first.name, f.name));
    const fieldCount = fields.length;
    if (fieldCount === 0) return;

    const precisions = getFloatPrecisions(first.name, fields, formatter);

    writeLine(builder, 'loop_');
    for (let i = 0; i < fieldCount; i++) {
        writeLine(builder, `_${first.name}.${fields[i].name}`);
    }

    for (let _c = 0; _c < categories.length; _c++) {
        const category = categories[_c];
        const data = category.data;

        if (category.rowCount === 0) continue;

        const it = category.keys ? category.keys() : Iterator.Range(0, category.rowCount - 1);
        while (it.hasNext)  {
            const key = it.move();

            let multiline = false;
            for (let _f = 0; _f < fieldCount; _f++) {
                multiline = writeValue(builder, data, key, fields[_f], precisions[_f]);
            }
            if (!multiline) StringBuilder.newline(builder);
        }
    }
    StringBuilder.write(builder, '#\n');
}

function isMultiline(value: string) {
    return !!value && value.indexOf('\n') >= 0;
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

    let escape = false, escapeCharStart = '\'', escapeCharEnd = '\' ';
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
    if (!escape && (fst === 35 /* # */ || fst === 59 /* ; */ || hasWhitespace)) {
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
