/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Table } from 'mol-data/db'
import Iterator from 'mol-data/iterator'
import * as Encoder from 'mol-io/writer/cif'

function columnValue(k: string) {
    return (i: number, d: any) => d[k].value(i);
}

function columnValueKind(k: string) {
    return (i: number, d: any) => d[k].valueKind(i);
}

function ofSchema(schema: Table.Schema) {
    const fields: Encoder.FieldDefinition[] = [];
    for (const k of Object.keys(schema)) {
        const t = schema[k];
        const type: any = t.valueKind === 'str' ? Encoder.FieldType.Str : t.valueKind === 'int' ? Encoder.FieldType.Int : Encoder.FieldType.Float;
        fields.push({ name: k, type, value: columnValue(k), valueKind: columnValueKind(k) })
    }
    return fields;
}

function ofTable<S extends Table.Schema>(name: string, table: Table<S>): Encoder.CategoryDefinition<number> {
    return { name, fields: ofSchema(table._schema) }
}

export function getCategoryInstanceProvider(name: string, table: Table<any>): Encoder.CategoryProvider {
    return () => {
        return {
            data: table,
            definition: ofTable(name, table),
            keys: () => Iterator.Range(0, table._rowCount - 1),
            rowCount: table._rowCount
        };
    }
}
