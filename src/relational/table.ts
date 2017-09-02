
import { Column } from './column'

/**
 * Represents that CIF category with multiple fields represented as columns.
 *
 * Example:
 * _category.field1
 * _category.field2
 * ...
 */
export interface Table {
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
    getColumn(name: string): Column;
}
