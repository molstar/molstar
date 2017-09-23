/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnType, createArray } from '../../column'
import { trimStr } from '../tokenizer'
import { parseIntSkipLeadingWhitespace, parseFloatSkipLeadingWhitespace } from '../number-parser'

export interface FixedColumnInfo {
    data: string,
    lines: ArrayLike<number>,
    rowCount: number
}

export default function FixedColumnProvider(info: FixedColumnInfo) {
    return function<T extends ColumnType>(offset: number, width: number, type: T) {
        return FixedColumn(info, offset, width, type);
    }
}

function getArrayValues(value: (row: number) => any, target: any[], start: number) {
    for (let i = 0, _e = target.length; i < _e; i++) target[i] = value(start + i);
    return target;
}

export function FixedColumn<T extends ColumnType>(info: FixedColumnInfo, offset: number, width: number, type: T): Column<T['@type']> {
    const { data, lines, rowCount } = info;
    const { kind } = type;

    const value: Column<T['@type']>['value'] = kind === 'str' ? row => {
        let s = lines[2 * row] + offset, e = s + width, le = lines[2 * row + 1];
        if (s >= le) return '';
        if (e > le) e = le;
        return trimStr(data, s, e);
    } : kind === 'int' ? row => {
        const s = lines[2 * row] + offset, e = s + width;
        return parseIntSkipLeadingWhitespace(data, s, e);
    } : row => {
        const s = lines[2 * row] + offset, e = s + width;
        return parseFloatSkipLeadingWhitespace(data, s, e);
    }
    return {
        isColumnDefined: true,
        rowCount,
        value,
        toArray(ctor, s, e) {
            const { array, start } = createArray(rowCount, ctor, s, e);
            return getArrayValues(value, array, start);
        }
    };
}