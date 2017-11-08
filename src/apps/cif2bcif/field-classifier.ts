/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-data/db'
import { Field } from 'mol-io/reader/cif/data-model'
import { FieldDefinition, FieldType } from 'mol-io/writer/cif/encoder'

const intRegex = /^-?\d+$/
const floatRegex = /^-?(([0-9]+)[.]?|([0-9]*[.][0-9]+))([(][0-9]+[)])?([eE][+-]?[0-9]+)?/

function classify(name: string, field: Field): FieldDefinition {
    let floatCount = 0, hasString = false;
    for (let i = 0, _i = field.rowCount; i < _i; i++) {
        const k = field.valueKind(i);
        if (k !== Column.ValueKind.Present) continue;
        const v = field.str(i);
        if (intRegex.test(v)) continue;
        else if (floatRegex.test(v)) floatCount++;
        else { hasString = true; break; }
    }

    if (hasString) return { name, type: FieldType.Str, value: field.str, valueKind: field.valueKind };
    if (floatCount > 0) return { name, type: FieldType.Float, value: field.float, valueKind: field.valueKind };
    return { name, type: FieldType.Int, value: field.int, valueKind: field.valueKind };
}

export default classify;