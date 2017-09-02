
import { ValuePresence } from './constants'

/**
 * A columns represents a single field of a CIF category.
 */
export interface Column {
    isDefined: boolean;

    getString(row: number): string | null;
    getInteger(row: number): number;
    getFloat(row: number): number;

    getValuePresence(row: number): ValuePresence;

    areValuesEqual(rowA: number, rowB: number): boolean;
    stringEquals(row: number, value: string): boolean;
}

/**
 * Represents a column that is not present.
 */
class _UndefinedColumn implements Column {  // tslint:disable-line:class-name
    isDefined = false;
    getString(row: number): string | null { return null; };
    getInteger(row: number): number { return 0; }
    getFloat(row: number): number { return 0.0; }
    getValuePresence(row: number): ValuePresence { return ValuePresence.NotSpecified; }
    areValuesEqual(rowA: number, rowB: number): boolean { return true; }
    stringEquals(row: number, value: string): boolean { return value === null; }
}
export const UndefinedColumn = new _UndefinedColumn() as Column;
