/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Table from './table'

/** A collection of tables */
type Database<Schema extends Database.Schema> = {
    readonly _name: string,
    readonly _tableNames: string[],
    readonly _schema: Schema
} & Database.Tables<Schema>

namespace Database {
    export type Tables<S extends Schema> = { [T in keyof Schema]: Table<S[T]> }
    export type Schema = { [table: string]: Table.Schema }

    export function ofTables<S extends Schema, Db = Database<S>>(name: string, schema: Schema, tables: Tables<S>): Db {
        const keys = Object.keys(tables);
        const ret = Object.create(null);
        const tableNames: string[] = [];
        ret._name = name;
        ret._tableNames = tableNames;
        ret._schema = schema;
        for (const k of keys) {
            if (!Table.is(tables[k])) continue;
            ret[k] = tables[k];
            tableNames[tableNames.length] = k;
        }
        return ret;
    }
}

export default Database