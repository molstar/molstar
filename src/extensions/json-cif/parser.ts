/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnHelpers } from '../../mol-data/db';
import { CifBlock, CifCategory, CifField, CifFile } from '../../mol-io/reader/cif';
import { ReaderResult } from '../../mol-io/reader/result';
import { Task } from '../../mol-task';
import { JSONCifCategory, JSONCifFile } from './model';

function Field(rows: Record<string, any>[], name: string): CifField {
    const str: CifField['str'] = row => {
        const v = rows[row][name];
        if (v === null || v === undefined) return '';
        if (typeof v === 'string') return v;
        return '' + v;
    };

    const number: CifField['int'] = row => +rows[row][name];
    const valueKind: CifField['valueKind'] = row => {
        const v = rows[row][name];
        if (v === null) return Column.ValueKinds.NotPresent;
        if (v === undefined) return Column.ValueKinds.Unknown;
        return Column.ValueKinds.Present;
    };

    const rowCount = rows.length;
    return {
        __array: undefined,
        binaryEncoding: undefined,
        isDefined: true,
        rowCount,
        str,
        int: number,
        float: number,
        valueKind,
        areValuesEqual: (rowA, rowB) => rows[rowA][name] === rows[rowB][name],
        toStringArray: params => ColumnHelpers.createAndFillArray(rowCount, str, params),
        toIntArray: params => ColumnHelpers.createAndFillArray(rowCount, number, params),
        toFloatArray: params => ColumnHelpers.createAndFillArray(rowCount, number, params),
    };
}

function Category(data: JSONCifCategory): CifCategory {
    const nameSet = new Set(data.fieldNames);
    const cache: Record<string, CifField> = Object.create(null);

    return {
        rowCount: data.rows.length,
        name: data.name,
        fieldNames: data.fieldNames,
        getField(name) {
            if (!nameSet.has(name)) return void 0;
            if (!!cache[name]) return cache[name];
            cache[name] = Field(data.rows, name);
            return cache[name];
        }
    };
}

function checkVersions(min: number[], current: number[]) {
    for (let i = 0; i < 2; i++) {
        if (min[i] > current[i]) return false;
    }
    return true;
}

export function parseJSONCif(data: JSONCifFile) {
    const minVersion = [0, 1];

    if (!checkVersions(minVersion, data.version.match(/(\d)\.(\d)\.\d/)!.slice(1).map(v => +v))) {
        throw new Error(`Unsupported format version. Current ${data.version}, required ${minVersion.join('.')}.`);
    }

    return CifFile(data.dataBlocks.map(block => {
        const cats = Object.create(null);
        for (const cat of block.categoryNames) cats[cat] = Category(block.categories[cat]);
        return CifBlock(block.categoryNames, cats, block.header);
    }));
}

export function parseJSONCifString(data: string) {
    return Task.create<ReaderResult<CifFile>>('Parse BinaryCIF', async ctx => {
        try {
            const json = JSON.parse(data) as JSONCifFile;
            const file = parseJSONCif(json);
            return ReaderResult.success(file);
        } catch (e) {
            return ReaderResult.error<CifFile>('' + e);
        }
    });
}