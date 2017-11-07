/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnHelpers } from 'mol-data/db'
import { Tokens } from '../tokenizer'
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '../number-parser'

export default function TokenColumnProvider(tokens: Tokens) {
    return function<T extends Column.Type>(type: T) {
        return TokenColumn(tokens, type);
    }
}

export function TokenColumn<T extends Column.Type>(tokens: Tokens, type: T): Column<T['T']> {
    const { data, indices, count: rowCount } = tokens;
    const { kind } = type;

    const value: Column<T['T']>['value'] =
          kind === 'str'
        ? row => data.substring(indices[2 * row], indices[2 * row + 1])
        : kind === 'int'
        ? row => fastParseInt(data, indices[2 * row], indices[2 * row + 1]) || 0
        : row => fastParseFloat(data, indices[2 * row], indices[2 * row + 1]) || 0;

    return {
        '@type': type,
        '@array': void 0,
        isDefined: true,
        rowCount,
        value,
        valueKind: row => Column.ValueKind.Present,
        toArray: params => ColumnHelpers.createAndFillArray(rowCount, value, params),
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