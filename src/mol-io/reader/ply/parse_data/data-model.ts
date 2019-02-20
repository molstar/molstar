/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifField as PlyColumn } from '../../cif/data-model'

export { PlyColumn }

export interface PlyFile {
    readonly name?: string,
    readonly PLY_File: ply_form
}

export function PlyFile(PLY_File: ply_form, name?: string): PlyFile {
    return { name, PLY_File };
}

export interface ply_form {
    readonly vertexCount: number,
    readonly faceCount: number,
    readonly propertyCount: number,
    readonly initialHead: ReadonlyArray<string>,
    readonly propertyNames: ReadonlyArray<string>,
    readonly properties: number[],
    readonly vertices: number[],
    readonly colors: number[],
    readonly normals: number[],
    readonly  faces: number[],
}

export function PlyStructure(vertexCount: number, faceCount: number, propertyCount: number, initialHead: string[], propertyNames: string[],  properties: number[],  vertices: number[],  colors: number[],  normals: number[], faces: number[]): ply_form {
    return {vertexCount, faceCount, propertyCount, initialHead: [...initialHead], propertyNames: [...propertyNames], properties: [...properties], vertices: [...vertices], colors: [...colors], normals: [...normals], faces: [...faces]};
}

export type PlyColumns = { [name: string]: PlyColumn }

// export namespace CsvTable {
//     export function empty(name: string): Table {
//         return { rowCount: 0, name, fieldNames: [], getColumn(name: string) { return void 0; } };
//     };
// }