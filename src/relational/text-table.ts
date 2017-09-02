

import { Table } from './table'
import { UndefinedColumn } from './column'
import { TextColumn, CifColumn } from './text-column'

import { Tokens } from '../utils/tokens'

/**
 * Represents a table backed by a string.
 */
export class TextTable implements Table {
    protected data: string;
    protected columnNameList: string[];
    protected columnIndices: Map<string, number>;

    /**
     * Name of the category.
     */
    name: string;

    /**
     * The array of columns.
     */
    get columnNames() {
        return this.columnNameList;
    }

    /**
     * Number of columns in the category.
     */
    columnCount: number;

    /**
     * Number of rows in the category.
     */
    rowCount: number;

    /**
     * Pairs of (start at index 2 * i, end at index 2 * i + 1) indices to the data string.
     * The "end" character is not included (for it's iterated as for (i = start; i < end; i++)).
     */
    indices: Int32Array;

    /**
     * Get a column object that makes accessing data easier.
     */
    getColumn(name: string): TextColumn {
        let i = this.columnIndices.get(name);
        if (i !== void 0) return new TextColumn(this, this.data, name, i);
        return UndefinedColumn as TextColumn;
    }

    initColumns(columns: string[]): void {
        this.columnIndices = new Map<string, number>();
        this.columnNameList = [];
        for (let i = 0; i < columns.length; i++) {
            this.columnIndices.set(columns[i], i);
            this.columnNameList.push(columns[i]);
        }
    }

    constructor(
        data: string, name: string, columns: string[], tokens: Tokens) {
        this.name = name;
        this.indices = tokens.indices;
        this.data = data;

        this.columnCount = columns.length;
        this.rowCount = (tokens.count / 2 / columns.length) | 0;

        this.initColumns(columns)
    }
}

export class CifTable extends TextTable {
    getColumn(name: string): CifColumn {
        let i = this.columnIndices.get(name);
        if (i !== void 0) return new CifColumn(this, this.data, name, i);
        return UndefinedColumn as CifColumn;
    }

    initColumns(columns: string[]): void {
        this.columnIndices = new Map<string, number>();
        this.columnNameList = [];
        for (let i = 0; i < columns.length; i++) {
            let colName = columns[i].substr(this.name.length + 1);
            this.columnIndices.set(colName, i);
            this.columnNameList.push(colName);
        }
    }
}
