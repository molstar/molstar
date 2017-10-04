/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnType, createAndFillArray } from '../../column'
import { Tokens } from '../tokenizer'
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '../number-parser'
import StringPool from '../../../../utils/short-string-pool'

export default function TokenColumnProvider(tokens: Tokens) {
    return function<T extends ColumnType>(type: T) {
        return TokenColumn(tokens, type);
    }
}

export function TokenColumn<T extends ColumnType>(tokens: Tokens, type: T): Column<T['@type']> {
    const { data, indices, count: rowCount } = tokens;
    const { kind } = type;
    const pool = kind === 'pooled-str' ? StringPool.create() : void 0;

    const value: Column<T['@type']>['value'] =
          kind === 'str'
        ? row => data.substring(indices[2 * row], indices[2 * row + 1])
        : kind === 'pooled-str'
        ? row => StringPool.get(pool!, data.substring(indices[2 * row], indices[2 * row + 1]))
        : kind === 'int'
        ? row => fastParseInt(data, indices[2 * row], indices[2 * row + 1]) || 0
        : row => fastParseFloat(data, indices[2 * row], indices[2 * row + 1]) || 0;

    return {
        isDefined: true,
        rowCount,
        value,
        isValueDefined(row) { return true; },
        toArray(params) { return createAndFillArray(rowCount, value, params); },
        stringEquals(row, v) {
            const s = indices[2 * row];
            const value = v || '';
            const len = value.length;
            if (len !== indices[2 * row + 1] - s) return false;
            for (let i = 0; i < len; i++) {
                if (data.charCodeAt(i + s) !== value.charCodeAt(i)) return false;
            }
            return true;
        },
        areValuesEqual: areValuesEqualProvider(tokens)
    };
}

export function areValuesEqualProvider(tokens: Tokens) {
    const { data, indices } = tokens;
    return function(rowA: number, rowB: number) {
        const aS = indices[2 * rowA], bS = indices[2 * rowB];
        const len = indices[2 * rowA + 1] - aS;
        if (len !== indices[2 *  rowB + 1] - bS) return false;
        for (let i = 0; i < len; i++) {
            if (data.charCodeAt(i + aS) !== data.charCodeAt(i + bS)) {
                return false;
            }
        }
        return true;
    }
}