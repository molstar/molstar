/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ReaderResult as Result } from '../result';
import { Task } from '../../../mol-task';
import { parseCsv } from '../csv/parser';
import { Column, Table } from '../../../mol-data/db';
import { toTable } from '../cif/schema';

import Schema = Column.Schema
import { CsvTable } from '../csv/data-model';


export const Schema3DG = {
    /** Chromosome name */
    chromosome: Schema.str,
    /** Base position */
    position: Schema.int,
    /** X coordinate */
    x: Schema.float,
    /** Y coordinate */
    y: Schema.float,
    /** Z coordinate */
    z: Schema.float,
};
export type Schema3DG = typeof Schema3DG

export interface File3DG {
    table: Table<Schema3DG>
}

const FieldNames = [ 'chromosome', 'position', 'x', 'y', 'z' ];

function categoryFromTable(name: string, table: CsvTable) {
    return {
        name,
        rowCount: table.rowCount,
        fieldNames: FieldNames,
        getField: (name: string) => {
            return table.getColumn(FieldNames.indexOf(name).toString());
        }
    };
}

export function parse3DG(data: string) {
    return Task.create<Result<File3DG>>('Parse 3DG', async ctx => {
        const opts = { quote: '', comment: '#', delimiter: '\t', noColumnNames: true };
        const csvFile = await parseCsv(data, opts).runInContext(ctx);
        if (csvFile.isError) return Result.error(csvFile.message, csvFile.line);
        const category = categoryFromTable('3dg', csvFile.result.table);
        const table = toTable(Schema3DG, category);
        return Result.success({ table });
    });
}