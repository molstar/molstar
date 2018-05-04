/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifField as CsvColumn } from '../cif/data-model'

export { CsvColumn }

export interface CsvFile {
    readonly name?: string,
    readonly table: CsvTable
}

export function CsvFile(table: CsvTable, name?: string): CsvFile {
    return { name, table };
}

export interface CsvTable {
    readonly rowCount: number,
    readonly columnNames: ReadonlyArray<string>,
    getColumn(name: string): CsvColumn | undefined
}

export function CsvTable(rowCount: number, columnNames: string[], columns: CsvColumns): CsvTable {
    return { rowCount, columnNames: [...columnNames], getColumn(name) { return columns[name]; } };
}

export type CsvColumns = { [name: string]: CsvColumn }

// export namespace CsvTable {
//     export function empty(name: string): Table {
//         return { rowCount: 0, name, fieldNames: [], getColumn(name: string) { return void 0; } };
//     };
// }