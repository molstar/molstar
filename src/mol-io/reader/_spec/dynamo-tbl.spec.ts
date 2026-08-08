/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author OpenAI
 */

import { Column } from '../../../mol-data/db';
import { parseDynamoTbl } from '../dynamo/tbl';

test('parses Dynamo TBL rows into named typed columns', async () => {
    const data = [
        '1 0 0 10 20 30 90 0 0 0 0 0 0 0 0 0 0 0 0 7 8 9 10 100 200 300 0 0 0 0 0 0 0 0 0 6.5',
        '2 0 0 1 2 3 0 90 0 0 0 0 0 0 0 0 0 0 0 7 8 9 10 40 50 60 0 0 0 0 0 0 0 0 0 6.5',
    ].join('\n');

    const parsed = await parseDynamoTbl(data).run();
    if (parsed.isError) throw new Error(parsed.message);

    const { fields, rowCount } = parsed.result;
    expect(rowCount).toBe(2);
    expect(fields.x.value(0)).toBe(100);
    expect(fields.y.value(0)).toBe(200);
    expect(fields.z.value(0)).toBe(300);
    expect(fields.tdrot.value(0)).toBe(90);
    expect(fields.tilt.value(1)).toBe(90);
    expect(fields.tomo.value(0)).toBe(7);
    expect(fields.apix.value(0)).toBe(6.5);
    expect(parsed.warnings).toEqual([]);
});

test('accepts Dynamo TBL rows with the 26 required columns', async () => {
    const tokens = ['1', '0', '0', '10', '20', '30', '45', '30', '60'];
    while (tokens.length < 26) tokens.push('0');
    tokens[23] = '100';
    tokens[24] = '200';
    tokens[25] = '300';
    const parsed = await parseDynamoTbl(tokens.join(' ')).run();
    if (parsed.isError) throw new Error(parsed.message);

    const { fields, rowCount } = parsed.result;
    expect(rowCount).toBe(1);
    expect(fields.tag.value(0)).toBe(1);
    expect(fields.dx.value(0)).toBe(10);
    expect(fields.narot.value(0)).toBe(60);
    expect(fields.x.value(0)).toBe(100);
    expect(fields.z.value(0)).toBe(300);
    // Optional columns beyond the 26 required ones are absent.
    expect(fields.apix.valueKind(0)).toBe(Column.ValueKinds.NotPresent);
    expect(fields.eig1.valueKind(0)).toBe(Column.ValueKinds.NotPresent);
});

test('rejects Dynamo TBL rows with fewer than 26 columns', async () => {
    const parsed = await parseDynamoTbl('1 0 0 10 20 30 45 30 60').run();
    expect(parsed.isError).toBe(true);
});
