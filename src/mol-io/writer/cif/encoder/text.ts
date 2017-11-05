/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from CIFTools.js (https://github.com/dsehnal/CIFTools.js)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-base/collections/database'
import StringBuilder from 'mol-base/utils/string-builder'
import * as Enc from '../encoder'
import Writer from '../../writer'

export default class TextCIFEncoder<Context> implements Enc.CIFEncoder<string, Context> {
    private builder = StringBuilder.create();
    private encoded = false;
    private dataBlockCreated = false;

    startDataBlock(header: string) {
        this.dataBlockCreated = true;
        StringBuilder.write(this.builder, `data_${(header || '').replace(/[ \n\t]/g, '').toUpperCase()}\n#\n`);
    }

    writeCategory(category: Enc.CategoryProvider, contexts?: Context[]) {
        if (this.encoded) {
            throw new Error('The writer contents have already been encoded, no more writing.');
        }

        if (!this.dataBlockCreated) {
            throw new Error('No data block created.');
        }

        const categories = !contexts || !contexts.length ? [category(<any>void 0)] : contexts.map(c => category(c));
        if (!categories.length) return;

        const rowCount = categories.reduce((v, c) => v + c.rowCount, 0);

        if (rowCount === 0) return;

        if (rowCount === 1) {
            writeCifSingleRecord(categories[0]!, this.builder);
        } else {
            writeCifLoop(categories, this.builder);
        }
    }

    encode() {
        this.encoded = true;
    }

    writeTo(stream: Writer<string>) {
        const chunks = StringBuilder.getChunks(this.builder);
        for (let i = 0, _i = chunks.length; i < _i; i++) {
            stream.write(chunks[i]);
        }
    }

    getData() {
        return StringBuilder.getString(this.builder);
    }
}

function writeValue(builder: StringBuilder, data: any, key: any, f: Enc.FieldDefinition<any, any>): boolean {
    const kind = f.valueKind;
    const p = kind ? kind(key, data) : Column.ValueKind.Present;
    if (p !== Column.ValueKind.Present) {
        if (p === Column.ValueKind.NotPresent) writeNotPresent(builder);
        else writeUnknown(builder);
    } else {
        const val = f.value(key, data);
        const t = f.type;
        if (t === Enc.FieldType.Str) {
            if (isMultiline(val as string)) {
                writeMultiline(builder, val as string);
                return true;
            } else {
                return writeChecked(builder, val as string);
            }
        } else if (t === Enc.FieldType.Int) {
            writeInteger(builder, val as number);
        } else {
            writeFloat(builder, val as number, 1000000);
        }
    }
    return false;
}

function writeCifSingleRecord(category: Enc.CategoryInstance<any>, builder: StringBuilder) {
    const fields = category.definition.fields;
    const data = category.data;
    const width = fields.reduce((w, s) => Math.max(w, s.name.length), 0) + category.definition.name.length + 6;

    const it = category.keys();
    const key = it.move();

    for (let _f = 0; _f < fields.length; _f++) {
        const f = fields[_f];
        StringBuilder.writePadRight(builder, `_${category.definition.name}.${f.name}`, width);
        const multiline = writeValue(builder, data, key, f);
        if (!multiline) StringBuilder.newline(builder);
    }
    StringBuilder.write(builder, '#\n');
}

function writeCifLoop(categories: Enc.CategoryInstance[], builder: StringBuilder) {
    const first = categories[0];
    const fields = first.definition.fields;
    const fieldCount = fields.length;

    writeLine(builder, 'loop_');
    for (let i = 0; i < fieldCount; i++) {
        writeLine(builder, `_${first.definition.name}.${fields[i].name}`);
    }

    for (let _c = 0; _c < categories.length; _c++) {
        const category = categories[_c];
        const data = category.data;

        if (category.rowCount === 0) continue;

        const it = category.keys();
        while (it.hasNext)  {
            const key = it.move();

            let multiline = false;
            for (let _f = 0; _f < fieldCount; _f++) {
                multiline = writeValue(builder, data, key, fields[_f]);
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
