/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from './column'
import { ValuePresence } from './constants'
import { TextTable } from './text-table'

import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '../utils/number-parser'
import { ShortStringPool } from '../utils/short-string-pool'

/**
 * Represents a single column.
 */
export class TextColumn implements Column {

    protected indices: Int32Array;
    protected columnCount: number;
    protected rowCount: number;
    protected stringPool = ShortStringPool.create();

    isDefined = true;

    /**
     * Returns the string value at given row.
     */
    getString(row: number): string | null {
        let i = (row * this.columnCount + this.index) * 2;
        return ShortStringPool.get(this.stringPool, this.data.substring(this.indices[i], this.indices[i + 1]));
    }

    /**
     * Returns the integer value at given row.
     */
    getInteger(row: number): number {
        let i = (row * this.columnCount + this.index) * 2;
        return fastParseInt(this.data, this.indices[i], this.indices[i + 1]);
    }

    /**
     * Returns the float value at given row.
     */
    getFloat(row: number): number {
        let i = (row * this.columnCount + this.index) * 2;
        return fastParseFloat(this.data, this.indices[i], this.indices[i + 1]);
    }

    /**
     * Returns true if the token has the specified string value.
     */
    stringEquals(row: number, value: string) {
        let aIndex = (row * this.columnCount + this.index) * 2,
            s = this.indices[aIndex],
            len = value.length;
        if (len !== this.indices[aIndex + 1] - s) return false;
        for (let i = 0; i < len; i++) {
            if (this.data.charCodeAt(i + s) !== value.charCodeAt(i)) return false;
        }
        return true;
    }

    /**
     * Determines if values at the given rows are equal.
     */
    areValuesEqual(rowA: number, rowB: number): boolean {
        const aIndex = (rowA * this.columnCount + this.index) * 2
        const bIndex = (rowB * this.columnCount + this.index) * 2
        const aS = this.indices[aIndex]
        const bS = this.indices[bIndex]
        const len = this.indices[aIndex + 1] - aS
        if (len !== this.indices[bIndex + 1] - bS) return false;
        for (let i = 0; i < len; i++) {
            if (this.data.charCodeAt(i + aS) !== this.data.charCodeAt(i + bS)) {
                return false;
            }
        }
        return true;
    }

    getValuePresence(row: number): ValuePresence {
        let index = 2 * (row * this.columnCount + this.index);
        if (this.indices[index] === this.indices[index + 1]) {
            return ValuePresence.NotSpecified
        }
        return ValuePresence.Present
    }

    constructor(table: TextTable, protected data: string, public name: string, public index: number) {
        this.indices = table.indices;
        this.columnCount = table.columnCount;
    }
}

export class CifColumn extends TextColumn {
    /**
     * Returns the string value at given row.
     */
    getString(row: number): string | null {
        let ret = super.getString(row)
        if (ret === '.' || ret === '?') return null;
        return ret;
    }

    /**
     * Returns true if the value is not defined (. or ? token).
     */
    getValuePresence(row: number): ValuePresence {
        let index = 2 * (row * this.columnCount + this.index);
        let s = this.indices[index];
        if (this.indices[index + 1] - s !== 1) return ValuePresence.Present;
        let v = this.data.charCodeAt(s);
        if (v === 46 /* . */) return ValuePresence.NotSpecified;
        if (v === 63 /* ? */) return ValuePresence.Unknown;
        return ValuePresence.Present;
    }
}
