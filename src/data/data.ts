/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export interface File {
    readonly name?: string,
    readonly blocks: ReadonlyArray<Block>
}

export interface Block {
    readonly header?: string,
    readonly categories: { readonly [name: string]: Category }
}

export interface Category {
    readonly rowCount: number,
    getField(name: string): Field | undefined
}

export namespace Category {
    export const Empty: Category = { rowCount: 0, getField(name: string) { return void 0; } };
}

export const enum ValuePresence {
    Present = 0,
    NotSpecified = 1,
    Unknown = 2
}

export type FieldArray = number[] | Float32Array | Float64Array | Int8Array | Int16Array | Int32Array | Uint8Array | Uint16Array | Uint32Array

/**
 * Implementation note:
 * Always implement this as a "plain" object so that the functions are "closures"
 * by default. This is to ensure that the schema access works without definiting
 * additional closures.
 */
export interface Field {
    readonly isDefined: boolean,

    str(row: number): string | null,
    int(row: number): number,
    float(row: number): number,
    bin(row: number): Uint8Array | null,

    presence(row: number): ValuePresence,

    areValuesEqual(rowA: number, rowB: number): boolean,
    stringEquals(row: number, value: string | null): boolean,

    toStringArray(startRow: number, endRowExclusive: number, ctor: (size: number) => FieldArray): ReadonlyArray<string>,
    toNumberArray(startRow: number, endRowExclusive: number, ctor: (size: number) => FieldArray): ReadonlyArray<number>
}