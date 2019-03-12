/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Schäfer, Marco <marco.schaefer@uni-tuebingen.de>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifField as PlyColumn } from '../cif/data-model'

export { PlyColumn }

export interface PlyFile {
    readonly name?: string,
    readonly PLY_File: PlyData
}

export function PlyFile(PLY_File: PlyData, name?: string): PlyFile {
    return { name, PLY_File };
}

export interface PlyData {
    readonly vertexCount: number,
    readonly faceCount: number,
    readonly propertyCount: number,
    readonly initialHead: ReadonlyArray<string>,
    readonly propertyNames: ReadonlyArray<string>,
    readonly properties: number[],
    readonly vertices: number[],
    readonly colors: number[],
    readonly normals: number[],
    readonly faces: number[],
}

// TODO note, removed `faces: [...faces]` pattern as that copies the data which I assume was not intentional (alex)
export function PlyData(vertexCount: number, faceCount: number, propertyCount: number, initialHead: string[], propertyNames: string[],  properties: number[],  vertices: number[],  colors: number[],  normals: number[], faces: number[]): PlyData {
    return {
        vertexCount,
        faceCount,
        propertyCount,
        initialHead,
        propertyNames,
        properties,
        vertices,
        colors,
        normals,
        faces
    };
}