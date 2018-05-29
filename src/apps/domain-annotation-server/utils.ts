/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Table } from 'mol-data/db'
import { CIFField, CIFCategory } from 'mol-io/writer/cif'

function columnValue(k: string) {
    return (i: number, d: any) => d[k].value(i);
}

function columnValueKind(k: string) {
    return (i: number, d: any) => d[k].valueKind(i);
}

function ofSchema(schema: Table.Schema) {
    const fields: CIFField[] = [];
    for (const k of Object.keys(schema)) {
        const t = schema[k];
        const type: any = t.valueType === 'str' ? CIFField.Type.Str : t.valueType === 'int' ? CIFField.Type.Int : CIFField.Type.Float;
        fields.push({ name: k, type, value: columnValue(k), valueKind: columnValueKind(k) })
    }
    return fields;
}

export function getCategoryInstanceProvider(name: string, table: Table<any>): CIFCategory.Provider {
    return () => {
        return {
            data: table,
            name,
            fields: ofSchema(table._schema),
            rowCount: table._rowCount
        };
    }
}
