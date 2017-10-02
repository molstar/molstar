/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnType, createArray } from '../../column'
import { trimStr, Lines } from '../tokenizer'
import { parseIntSkipLeadingWhitespace, parseFloatSkipLeadingWhitespace } from '../number-parser'
import StringPool from '../../../../utils/short-string-pool'

export default function FixedColumnProvider(lines: Lines) {
    return function<T extends ColumnType>(offset: number, width: number, type: T) {
        return FixedColumn(lines, offset, width, type);
    }
}

function fillArrayValues(value: (row: number) => any, target: any[], start: number) {
    for (let i = 0, _e = target.length; i < _e; i++) target[i] = value(start + i);
    return target;
}

export function FixedColumn<T extends ColumnType>(lines: Lines, offset: number, width: number, type: T): Column<T['@type']> {
    const { data, tokens, count: rowCount } = lines;
    const { kind } = type;
    const pool = kind === 'pooled-str' ? StringPool.create() : void 0;

    const value: Column<T['@type']>['value'] = kind === 'str' ? row => {
        let s = tokens[2 * row] + offset, le = tokens[2 * row + 1];
        if (s >= le) return '';
        let e = s + width;
        if (e > le) e = le;
        return trimStr(data, s, e);
    } : kind === 'pooled-str' ? row => {
        let s = tokens[2 * row] + offset, le = tokens[2 * row + 1];
        if (s >= le) return '';
        let e = s + width;
        if (e > le) e = le;
        return StringPool.get(pool!, trimStr(data, s, e));
    } : kind === 'int' ? row => {
        const s = tokens[2 * row] + offset;
        if (s > tokens[2 * row + 1]) return 0;
        return parseIntSkipLeadingWhitespace(data, s, s + width);
    } : row => {
        const s = tokens[2 * row] + offset;
        if (s > tokens[2 * row + 1]) return 0;
        return parseFloatSkipLeadingWhitespace(data, s, s + width);
    };
    return {
        isDefined: true,
        rowCount,
        value,
        isValueDefined(row) { return true; },
        toArray(params) {
            const { array, start } = createArray(rowCount, params);
            return fillArrayValues(value, array, start);
        },
        areValuesEqual(rowA, rowB) {
            return value(rowA) === value(rowB);
        }
    };
}