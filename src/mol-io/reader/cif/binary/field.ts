/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, ColumnHelpers } from '../../../../mol-data/db';
import * as Data from '../data-model';
import { EncodedColumn, decode } from '../../../common/binary-cif';
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '../../common/text/number-parser';

export default function Field(column: EncodedColumn): Data.CifField {
    const mask = column.mask ? decode(column.mask) as number[] : void 0;
    const data = decode(column.data);
    const isNumeric = ColumnHelpers.isTypedArray(data);

    const str: Data.CifField['str'] = isNumeric
        ? mask
            ? row => mask[row] === Column.ValueKind.Present ? '' + data[row] : ''
            : row => '' + data[row]
        : mask
            ? row => mask[row] === Column.ValueKind.Present ? data[row] : ''
            : row => data[row];

    const int: Data.CifField['int'] = isNumeric
        ? row => data[row]
        : row => { const v = data[row]; return fastParseInt(v, 0, v.length); };

    const float: Data.CifField['float'] = isNumeric
        ? row => data[row]
        : row => { const v = data[row]; return fastParseFloat(v, 0, v.length); };

    const valueKind: Data.CifField['valueKind'] = mask
        ? row => mask[row]
        : row => Column.ValueKind.Present;

    const rowCount = data.length;

    return {
        __array: data,
        binaryEncoding: column.data.encoding,
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
    };
}