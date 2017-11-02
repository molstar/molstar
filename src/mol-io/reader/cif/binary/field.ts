/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, ColumnHelpers } from 'mol-base/collections/table'
import * as Data from '../data-model'
import { EncodedColumn } from './encoding'
import decode from './decoder'
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '../../common/text/number-parser'

function wrap(o: Data.Field) {
    return { ...o };
}

export default function Field(column: EncodedColumn): Data.Field {
    const mask = column.mask ? decode(column.mask) as number[] : void 0;
    const data = decode(column.data);
    const isNumeric = ColumnHelpers.isTypedArray(data);

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

    return wrap({
        '@array': data,
        isDefined: true,
        rowCount,
        str,
        int,
        float,
        valueKind,
        areValuesEqual: (rowA, rowB) => data[rowA] === data[rowB],
        toStringArray: params => ColumnHelpers.createAndFillArray(rowCount, str, params),
        toIntArray: isNumeric
            ? params => ColumnHelpers.typedArrayWindow(data, params)
            : params => ColumnHelpers.createAndFillArray(rowCount, int, params),
        toFloatArray: isNumeric
            ? params => ColumnHelpers.typedArrayWindow(data, params)
            : params => ColumnHelpers.createAndFillArray(rowCount, float, params)
    });
}

// return wrap({
//     '@array': data,
//     isDefined: true,
//     rowCount,
//     str: str as any,
//     int,
//     float,
//     valueKind,
//     areValuesEqual: (rowA, rowB) => data[rowA] === data[rowB],
//     toStringArray: params => ColumnHelpers.createAndFillArray(rowCount, str, params),
//     toIntArray: isNumeric
//         ? params => ColumnHelpers.typedArrayWindow(data, params)
//         : params => ColumnHelpers.createAndFillArray(rowCount, int, params),
//     toFloatArray: isNumeric
//         ? params => ColumnHelpers.typedArrayWindow(data, params)
//         : params => ColumnHelpers.createAndFillArray(rowCount, float, params)
// });