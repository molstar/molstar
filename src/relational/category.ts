/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, UndefinedColumn } from './column'

/**
 * Represents a tabular category with multiple fields represented as columns.
 *
 * Example:
 * _category.field1
 * _category.field2
 * ...
 */
export abstract class Category {
    name: string;
    rowCount: number;
    columnCount: number;
    columnNames: string[];

    /**
     * If a field with the given name is not present, returns UndefinedColumn.
     *
     * Columns are accessed by their field name only, i.e.
     * _category.field is accessed by
     * category.getColumn('field')
     *
     * Note that columns are created on demand and there is some computational
     * cost when creating a new column. Therefore, if you need to reuse a column,
     * it is a good idea to cache it.
     */
    abstract getColumn(name: string): Column;

    getColumnsFromSchema<T extends object> (schema: T) {
        return CategoryColumns(this, schema)
    }
}

/**
 * Represents a category that is not present.
 */
class _UndefinedCategory extends Category {  // tslint:disable-line:class-name
    name: ''
    rowCount = 0
    columnCount = 0
    columnNames = []
    getColumn(name: string) { return UndefinedColumn }
}
export const UndefinedCategory = new _UndefinedCategory() as Category;


export type CategoryColumns<Columns extends string> = { readonly [name in Columns]: Column }
export function CategoryColumns<T extends object>(category: Category | undefined, columns: T): CategoryColumns<keyof T> {
    const ret = Object.create(null);
    if (!category) for (const c of Object.keys(columns)) ret[c] = UndefinedColumn;
    else for (const c of Object.keys(columns)) ret[c] = category.getColumn(c);
    return ret;
}
