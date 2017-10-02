/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Column from '../common/column'
import * as Data from './data-model'
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '../common/text/number-parser'
import StringPool from '../../utils/short-string-pool'

export default function CifTextField(data: string, tokens: ArrayLike<number>, rowCount: number): Data.Field {
    const stringPool = StringPool.create();

    const str: Data.Field['str'] = row => {
        const ret = StringPool.get(stringPool, data.substring(tokens[2 * row], tokens[2 * row + 1]));
        if (ret === '.' || ret === '?') return '';
        return ret;
    };

    const int: Data.Field['int'] = row => {
        return fastParseInt(data, tokens[2 * row], tokens[2 * row + 1]) || 0;
    };

    const float: Data.Field['float'] = row => {
        return fastParseFloat(data, tokens[2 * row], tokens[2 * row + 1]) || 0;
    };

    const presence: Data.Field['presence'] = row => {
        const s = tokens[2 * row];
        if (tokens[2 * row + 1] - s !== 1) return Data.ValuePresence.Present;
        const v = data.charCodeAt(s);
        if (v === 46 /* . */) return Data.ValuePresence.NotSpecified;
        if (v === 63 /* ? */) return Data.ValuePresence.Unknown;
        return Data.ValuePresence.Present;
    };

    return {
        isDefined: true,
        rowCount,
        str,
        int,
        float,
        presence,
        areValuesEqual(rowA, rowB) {
            const aS = tokens[2 * rowA], bS = tokens[2 * rowB];
            const len = tokens[2 * rowA + 1] - aS;
            if (len !== tokens[2 *  rowB + 1] - bS) return false;
            for (let i = 0; i < len; i++) {
                if (data.charCodeAt(i + aS) !== data.charCodeAt(i + bS)) {
                    return false;
                }
            }
            return true;
        },
        stringEquals(row, value) {
            const s = tokens[2 * row];
            if (!value) return presence(row) !== Data.ValuePresence.Present;
            const len = value.length;
            if (len !== tokens[2 * row + 1] - s) return false;
            for (let i = 0; i < len; i++) {
                if (data.charCodeAt(i + s) !== value.charCodeAt(i)) return false;
            }
            return true;
        },
        toStringArray(params) {
            const { array, start } = Column.createArray(rowCount, params);
            return fillArrayValues(str, array, start);
        },
        toIntArray(params) {
            const { array, start } = Column.createArray(rowCount, params);
            return fillArrayValues(int, array, start);
        },
        toFloatArray(params) {
            const { array, start } = Column.createArray(rowCount, params);
            return fillArrayValues(float, array, start);
        }
    }
}

function fillArrayValues(value: (row: number) => any, target: any[], start: number) {
    for (let i = 0, _e = target.length; i < _e; i++) target[i] = value(start + i);
    return target;
}