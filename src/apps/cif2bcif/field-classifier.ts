/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-data/db'
import { CifField } from 'mol-io/reader/cif/data-model'
import { CifWriter } from 'mol-io/writer/cif'
import { classifyFloatArray, classifyIntArray } from 'mol-io/common/binary-cif/classifier';

const intRegex = /^-?\d+$/
const floatRegex = /^-?(([0-9]+)[.]?|([0-9]*[.][0-9]+))([(][0-9]+[)])?([eE][+-]?[0-9]+)?$/

// Classify a cif field as str, int or float based the data it contains.
// To classify a field as int or float all items are checked.
function classify(name: string, field: CifField): CifWriter.Field {
    let floatCount = 0, hasString = false;
    for (let i = 0, _i = field.rowCount; i < _i; i++) {
        const k = field.valueKind(i);
        if (k !== Column.ValueKind.Present) continue;
        const v = field.str(i);
        if (intRegex.test(v)) continue;
        else if (floatRegex.test(v)) floatCount++;
        else { hasString = true; break; }
    }

    if (hasString) return { name, type: CifWriter.Field.Type.Str, value: field.str, valueKind: field.valueKind };
    if (floatCount > 0) {
        const encoder = classifyFloatArray(field.toFloatArray({ array: Float64Array }));
        return CifWriter.Field.float(name, field.float, { valueKind: field.valueKind, encoder, typedArray: Float64Array });
    } else {
        const encoder = classifyIntArray(field.toIntArray({ array: Int32Array }));
        return CifWriter.Field.int(name, field.int, { valueKind: field.valueKind, encoder, typedArray: Int32Array });
    }
}

export default classify;