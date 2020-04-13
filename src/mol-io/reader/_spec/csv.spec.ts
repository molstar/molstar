/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseCsv } from '../csv/parser';

const csvStringBasic = `StrCol,IntCol,FloatCol
# comment
string1,-1,-0.34e3
string2,42,2.44`;

const csvStringAdvanced = `StrCol,"Int Col",FloatCol
 string1  \t , -1,  -0.34e3
    # comment
   " stri
ng2" ,42, 2.44 `;

const tabString = `StrCol\tIntCol\tFloatCol
string1\t-1\t-0.34e3
string2\t42\t2.44`;

describe('csv reader', () => {
    it('basic', async () => {
        const parsed = await parseCsv(csvStringBasic).run();
        if (parsed.isError) return;
        const csvFile = parsed.result;

        // csvFile.table.columnNames.forEach(name => {
        //     const col = csvFile.table.getColumn(name)
        //     if (col) console.log(name, col.toStringArray())
        // })

        const strCol = csvFile.table.getColumn('StrCol');
        if (strCol) expect(strCol.toStringArray()).toEqual(['string1', 'string2']);

        const intCol = csvFile.table.getColumn('IntCol');
        if (intCol) expect(intCol.toIntArray()).toEqual([-1, 42]);

        const floatCol = csvFile.table.getColumn('FloatCol');
        if (floatCol) expect(floatCol.toFloatArray()).toEqual([-340.0, 2.44]);

        expect.assertions(3);
    });

    it('advanced', async () => {
        const parsed = await parseCsv(csvStringAdvanced).run();
        if (parsed.isError) return;
        const csvFile = parsed.result;

        const strCol = csvFile.table.getColumn('StrCol');
        if (strCol) expect(strCol.toStringArray()).toEqual(['string1', ' stri\nng2']);

        const intCol = csvFile.table.getColumn('Int Col');
        if (intCol) expect(intCol.toIntArray()).toEqual([-1, 42]);

        const floatCol = csvFile.table.getColumn('FloatCol');
        if (floatCol) expect(floatCol.toFloatArray()).toEqual([-340.0, 2.44]);

        expect.assertions(3);
    });

    it('tabs', async () => {
        const parsed = await parseCsv(tabString, { delimiter: '\t' }).run();
        if (parsed.isError) return;
        const csvFile = parsed.result;

        const strCol = csvFile.table.getColumn('StrCol');
        if (strCol) expect(strCol.toStringArray()).toEqual(['string1', 'string2']);

        const intCol = csvFile.table.getColumn('IntCol');
        if (intCol) expect(intCol.toIntArray()).toEqual([-1, 42]);

        const floatCol = csvFile.table.getColumn('FloatCol');
        if (floatCol) expect(floatCol.toFloatArray()).toEqual([-340.0, 2.44]);

        expect.assertions(3);
    });
});