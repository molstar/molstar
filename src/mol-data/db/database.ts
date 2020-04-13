/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Table from './table';

/** A collection of tables */
type Database<Schema extends Database.Schema> = {
    readonly _name: string,
    readonly _tableNames: ReadonlyArray<string>,
    readonly _schema: Schema
} & Database.Tables<Schema>

namespace Database {
    export type Tables<S extends Schema> = { [T in keyof S]: Table<S[T]> }
    export type Schema = { [table: string]: Table.Schema }

    export function ofTables<S extends Schema>(name: string, schema: Schema, tables: Tables<S>) {
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

    export function getTablesAsRows<S extends Schema>(database: Database<S>) {
        const ret: { [k: string]: Table.Row<any>[] } = {};
        for (const k of database._tableNames) {
            ret[k] = Table.getRows(database[k]);
        }
        return ret;
    }
}

export default Database;