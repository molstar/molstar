/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Field as Column } from '../cif/data-model'

export { Column }

export interface File {
    readonly name?: string,
    readonly table: Table
}

export function File(table: Table, name?: string): File {
    return { name, table };
}

export interface Table {
    readonly rowCount: number,
    readonly columnNames: ReadonlyArray<string>,
    getColumn(name: string): Column | undefined
}

export function Table(rowCount: number, columnNames: string[], columns: Columns): Table {
    return { rowCount, columnNames: [...columnNames], getColumn(name) { return columns[name]; } };
}

export type Columns = { [name: string]: Column }

// export namespace Table {
//     export function empty(name: string): Table {
//         return { rowCount: 0, name, fieldNames: [], getColumn(name: string) { return void 0; } };
//     };
// }