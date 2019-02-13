/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifField as CsvColumn } from '../../cif/data-model'

export { CsvColumn }

export interface PlyFile {
    readonly name?: string,
    readonly PLY_File: ply_form
}

export function CsvFile(PLY_File: ply_form, name?: string): PlyFile {
    return { name, PLY_File };
}

export interface ply_form {
    readonly rowCount: number,
    readonly vertexCount: number,
    readonly faceCount: number,
    readonly propertyCount: number,
    readonly initialHead: ReadonlyArray<string>,
    getColumn(name: string): CsvColumn | undefined
}

export function CsvTable(rowCount: number, vertexCount: number, faceCount: number, propertyCount: number, initialHead: string[], columns: CsvColumns): ply_form {
    return { rowCount, vertexCount, faceCount, propertyCount, initialHead: [...initialHead], getColumn(name) { return columns[name]; } };
}

export type CsvColumns = { [name: string]: CsvColumn }

// export namespace CsvTable {
//     export function empty(name: string): Table {
//         return { rowCount: 0, name, fieldNames: [], getColumn(name: string) { return void 0; } };
//     };
// }