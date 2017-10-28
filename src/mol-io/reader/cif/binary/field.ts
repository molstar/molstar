/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column, { isTypedArray, createAndFillArray, typedArrayWindow } from '../../../../mol-base/collections/column'
import * as Data from '../data-model'
import { EncodedColumn } from './encoding'
import decode from './decoder'
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '../../common/text/number-parser'

export default function Field(column: EncodedColumn): Data.Field {
    const mask = column.mask ? decode(column.mask) as number[] : void 0;
    const data = decode(column.data);
    const isNumeric = isTypedArray(data);

    const str: Data.Field['str'] = isNumeric
        ? mask
            ? row => mask[row] === Column.ValueKind.Present ? '' + data[row] : ''
            : row => '' + data[row]
        : mask
            ? row => mask[row] === Column.ValueKind.Present ? data[row] : ''
            : row => data[row];

    const int: Data.Field['int'] = isNumeric
        ? row => data[row]
        : row => { const v = data[row]; return fastParseInt(v, 0, v.length); };

    const float: Data.Field['float'] = isNumeric
        ? row => data[row]
        : row => { const v = data[row]; return fastParseFloat(v, 0, v.length); };

    const valueKind: Data.Field['valueKind'] = mask
        ? row => mask[row]
        : row => Column.ValueKind.Present;

    const rowCount = data.length;

    return {
        '@array': data,
        isDefined: true,
        rowCount,
        str,
        int,
        float,
        valueKind,
        areValuesEqual: (rowA, rowB) => data[rowA] === data[rowB],
        stringEquals: (row, v) => str(row) === v,
        toStringArray: params => createAndFillArray(rowCount, str, params),
        toIntArray: isNumeric
            ? params => typedArrayWindow(data, params)
            : params => createAndFillArray(rowCount, int, params),
        toFloatArray: isNumeric
            ? params => typedArrayWindow(data, params)
            : params => createAndFillArray(rowCount, float, params)
    };
}