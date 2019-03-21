/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from 'mol-data/db';

// http://paulbourke.net/dataformats/ply/
// https://en.wikipedia.org/wiki/PLY_(file_format)

export const PlyTypeByteLength = {
    'char': 1,
    'uchar': 1,
    'short': 2,
    'ushort': 2,
    'int': 4,
    'uint': 4,
    'float': 4,
    'double': 8
}
export type PlyType = keyof typeof PlyTypeByteLength
export const PlyTypes = new Set(Object.keys(PlyTypeByteLength))
export function PlyType(str: string) {
    if (!PlyTypes.has(str)) throw new Error(`unknown ply type '${str}'`)
    return str as PlyType
}

export interface PlyFile {
    readonly comments: ReadonlyArray<string>
    readonly elementNames: ReadonlyArray<string>
    getElement(name: string): PlyElement | undefined
}

export function PlyFile(elements: PlyElement[], elementNames: string[], comments: string[]): PlyFile {
    const elementMap = new Map<string, PlyElement>()
    for (let i = 0, il = elementNames.length; i < il; ++i) {
        elementMap.set(elementNames[i], elements[i])
    }
    return {
        comments,
        elementNames,
        getElement: (name: string) => {
            return elementMap.get(name)
        }
    };
}

export type PlyElement = PlyTable | PlyList

export interface PlyTable {
    readonly kind: 'table'
    readonly rowCount: number
    readonly propertyNames: ReadonlyArray<string>
    getProperty(name: string): Column<number> | undefined
}

export interface PlyListValue {
    readonly entries: ArrayLike<number>
    readonly count: number
}

export interface PlyList {
    readonly kind: 'list'
    readonly rowCount: number,
    readonly name: string,
    value: (row: number) => PlyListValue
}