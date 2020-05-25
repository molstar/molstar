/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface Database {
    tables: { [ tableName: string ]: Table }
    aliases: { [ path: string ]: string[] }
}
export interface Table {
    description: string
    key: Set<string>
    columns: { [ columnName: string ]: Column }
}
export type Column = IntCol | StrCol | FloatCol | CoordCol | EnumCol | VectorCol | MatrixCol | ListCol

type BaseCol = { description: string }

export type IntCol = { type: 'int' } & BaseCol
export function IntCol(description: string): IntCol { return { type: 'int', description }; }

export type StrCol = { type: 'str' } & BaseCol
export function StrCol(description: string): StrCol { return { type: 'str', description }; }

export type FloatCol = { type: 'float' } & BaseCol
export function FloatCol(description: string): FloatCol { return { type: 'float', description }; }

export type CoordCol = { type: 'coord' } & BaseCol
export function CoordCol(description: string): CoordCol { return { type: 'coord', description }; }

export type EnumCol = { type: 'enum', subType: 'int' | 'str', values: string[] } & BaseCol
export function EnumCol(values: string[], subType: 'int' | 'str', description: string): EnumCol {
    return { type: 'enum', description, values, subType };
}

export type VectorCol = { type: 'vector', length: number } & BaseCol
export function VectorCol(length: number, description: string): VectorCol {
    return { type: 'vector', description, length };
}

export type MatrixCol = { type: 'matrix', rows: number, columns: number } & BaseCol
export function MatrixCol(columns: number, rows: number, description: string): MatrixCol {
    return { type: 'matrix', description, columns, rows };
}

export type ListCol = { type: 'list', subType: 'int' | 'str' | 'float' | 'coord', separator: string } & BaseCol
export function ListCol(subType: 'int' | 'str' | 'float' | 'coord', separator: string, description: string): ListCol {
    return { type: 'list', description, separator, subType };
}

export type Filter = { [ table: string ]: { [ column: string ]: true } }

export function mergeFilters (...filters: Filter[]) {
    const n = filters.length;
    const mergedFilter: Filter = {};
    const fields: Map<string, number> = new Map();
    filters.forEach(filter => {
        Object.keys(filter).forEach(category => {
            Object.keys(filter[ category ]).forEach(field => {
                const key = `${category}.${field}`;
                const value = fields.get(key) || 0;
                fields.set(key, value + 1);
            });
        });
    });
    fields.forEach((v, k) => {
        if (v !== n) return;
        const [categoryName, fieldName] = k.split('.');
        if (categoryName in mergedFilter) {
            mergedFilter[categoryName][fieldName] = true;
        } else {
            mergedFilter[categoryName] = { fieldName: true };
        }
    });
    return mergedFilter;
}
