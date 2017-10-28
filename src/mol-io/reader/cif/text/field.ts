/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column, { createAndFillArray } from '../../../../mol-base/collections/column'
import * as TokenColumn from '../../common/text/column/token'
import { Tokens } from '../../common/text/tokenizer'
import * as Data from '../data-model'
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '../../common/text/number-parser'

export default function CifTextField(tokens: Tokens, rowCount: number): Data.Field {
    const { data, indices } = tokens;

    const str: Data.Field['str'] = row => {
        const ret = data.substring(indices[2 * row], indices[2 * row + 1]);
        if (ret === '.' || ret === '?') return '';
        return ret;
    };

    const int: Data.Field['int'] = row => {
        return fastParseInt(data, indices[2 * row], indices[2 * row + 1]) || 0;
    };

    const float: Data.Field['float'] = row => {
        return fastParseFloat(data, indices[2 * row], indices[2 * row + 1]) || 0;
    };

    const valueKind: Data.Field['valueKind'] = row => {
        const s = indices[2 * row];
        if (indices[2 * row + 1] - s !== 1) return Column.ValueKind.Present;
        const v = data.charCodeAt(s);
        if (v === 46 /* . */) return Column.ValueKind.NotPresent;
        if (v === 63 /* ? */) return Column.ValueKind.Unknown;
        return Column.ValueKind.Present;
    };

    return {
        '@array': void 0,
        isDefined: true,
        rowCount,
        str,
        int,
        float,
        valueKind,
        areValuesEqual: TokenColumn.areValuesEqualProvider(tokens),
        stringEquals: (row, v) => {
            const s = indices[2 * row];
            const value = v || '';
            if (!value && valueKind(row) !== Column.ValueKind.Present) return true;
            const len = value.length;
            if (len !== indices[2 * row + 1] - s) return false;
            for (let i = 0; i < len; i++) {
                if (data.charCodeAt(i + s) !== value.charCodeAt(i)) return false;
            }
            return true;
        },
        toStringArray: params => createAndFillArray(rowCount, str, params),
        toIntArray: params => createAndFillArray(rowCount, int, params),
        toFloatArray: params => createAndFillArray(rowCount, float, params)
    }
}