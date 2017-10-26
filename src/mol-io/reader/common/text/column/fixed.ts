/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnType, createAndFillArray } from '../../column'
import { trimStr, Tokens } from '../tokenizer'
import { parseIntSkipLeadingWhitespace, parseFloatSkipLeadingWhitespace } from '../number-parser'
import StringPool from '../../../../utils/short-string-pool'

export default function FixedColumnProvider(lines: Tokens) {
    return function<T extends ColumnType>(offset: number, width: number, type: T) {
        return FixedColumn(lines, offset, width, type);
    }
}

export function FixedColumn<T extends ColumnType>(lines: Tokens, offset: number, width: number, type: T): Column<T['@type']> {
    const { data, indices, count: rowCount } = lines;
    const { kind } = type;
    const pool = kind === 'pooled-str' ? StringPool.create() : void 0;

    const value: Column<T['@type']>['value'] = kind === 'str' ? row => {
        let s = indices[2 * row] + offset, le = indices[2 * row + 1];
        if (s >= le) return '';
        let e = s + width;
        if (e > le) e = le;
        return trimStr(data, s, e);
    } : kind === 'pooled-str' ? row => {
        let s = indices[2 * row] + offset, le = indices[2 * row + 1];
        if (s >= le) return '';
        let e = s + width;
        if (e > le) e = le;
        return StringPool.get(pool!, trimStr(data, s, e));
    } : kind === 'int' ? row => {
        const s = indices[2 * row] + offset;
        if (s > indices[2 * row + 1]) return 0;
        return parseIntSkipLeadingWhitespace(data, s, s + width);
    } : row => {
        const s = indices[2 * row] + offset;
        if (s > indices[2 * row + 1]) return 0;
        return parseFloatSkipLeadingWhitespace(data, s, s + width);
    };
    return {
        '@type': type,
        isDefined: true,
        rowCount,
        value,
        isValueDefined: row => true,
        toArray: params => createAndFillArray(rowCount, value, params),
        stringEquals: (row, v) => value(row) === v,
        areValuesEqual: (rowA, rowB) => value(rowA) === value(rowB)
    };
}