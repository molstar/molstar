/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from '../../../data/data'
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from './number-parser'
import { Tokens } from './tokenizer'
import ShortStringPool from '../../../utils/short-string-pool'

export function createCategory(data: string, fields: string[], tokens: Tokens, rowCount: number) {
    const fi: TokenFieldInfo = { data, fieldCount: fields.length, tokens: tokens.indices };

    const categoryFields = Object.create(null);
    for (let i = 0; i < fi.fieldCount; ++i) {
        categoryFields[fields[i]] = TokenField(fi, i);
    }
    return Data.Category(rowCount, categoryFields);
}

export interface TokenFieldInfo {
    data: string,
    tokens: ArrayLike<number>,
    fieldCount: number,
    isCif?: boolean
}

export function TokenField(info: TokenFieldInfo, index: number): Data.Field {
    const { data, tokens, fieldCount, isCif = false } = info;
    const stringPool = ShortStringPool.create();

    const str: Data.Field['str'] = isCif ? row => {
        const i = (row * fieldCount + index) * 2;
        const ret = ShortStringPool.get(stringPool, data.substring(tokens[i], tokens[i + 1]));
        if (ret === '.' || ret === '?') return null;
        return ret;
    } : row => {
        const i = (row * fieldCount + index) * 2;
        return ShortStringPool.get(stringPool, data.substring(tokens[i], tokens[i + 1]));
    };

    const int: Data.Field['int'] = row => {
        const i = (row * fieldCount + index) * 2;
        return fastParseInt(data, tokens[i], tokens[i + 1]) || 0;
    };

    const float: Data.Field['float'] = row => {
        const i = (row * fieldCount + index) * 2;
        return fastParseFloat(data, tokens[i], tokens[i + 1]) || 0;
    };

    const presence: Data.Field['presence'] = isCif ? row => {
        const i = 2 * (row * fieldCount + index);
        const s = tokens[i];
        if (tokens[i + 1] - s !== 1) return Data.ValuePresence.Present;
        const v = data.charCodeAt(s);
        if (v === 46 /* . */) return Data.ValuePresence.NotSpecified;
        if (v === 63 /* ? */) return Data.ValuePresence.Unknown;
        return Data.ValuePresence.Present;
    } : row => {
        const i = 2 * (row * fieldCount + index);
        return tokens[i] === tokens[i + 1] ? Data.ValuePresence.NotSpecified : Data.ValuePresence.Present
    };

    return {
        isDefined: true,
        str,
        int,
        float,
        value: str,
        presence,
        areValuesEqual: (rowA, rowB) => {
            const aI = (rowA * fieldCount + index) * 2, aS = tokens[aI];
            const bI = (rowB * fieldCount + index) * 2, bS = tokens[bI];
            const len = tokens[aI + 1] - aS;
            if (len !== tokens[bI + 1] - bS) return false;
            for (let i = 0; i < len; i++) {
                if (data.charCodeAt(i + aS) !== data.charCodeAt(i + bS)) {
                    return false;
                }
            }
            return true;
        },
        stringEquals: (row, value) => {
            const aI = (row * fieldCount + index) * 2;
            const s = tokens[aI];
            if (!value) return presence(row) !== Data.ValuePresence.Present;
            const len = value.length;
            if (len !== tokens[aI + 1] - s) return false;
            for (let i = 0; i < len; i++) {
                if (data.charCodeAt(i + s) !== value.charCodeAt(i)) return false;
            }
            return true;
        },
        toStringArray: (startRow, endRowExclusive, ctor) => {
            const count = endRowExclusive - startRow;
            const ret = ctor(count) as any;
            for (let i = 0; i < count; i++) { ret[i] = str(startRow + i); }
            return ret;
        },
        toIntArray: (startRow, endRowExclusive, ctor) => {
            const count = endRowExclusive - startRow;
            const ret = ctor(count) as any;
            for (let i = 0; i < count; i++) { ret[i] = int(startRow + i); }
            return ret;
        },
        toFloatArray: (startRow, endRowExclusive, ctor) => {
            const count = endRowExclusive - startRow;
            const ret = ctor(count) as any;
            for (let i = 0; i < count; i++) { ret[i] = float(startRow + i); }
            return ret;
        }
    }
}