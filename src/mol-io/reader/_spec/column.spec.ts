/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import FixedColumn from '../common/text/column/fixed';
import TokenColumn from '../common/text/column/token';
import { Column, ColumnHelpers } from '../../../mol-data/db';

const lines = [
    '1.123 abc',
    '1.00  a',
    '1.1    bcd   ',
    '',
    ' 5'
];

const linesData = lines.join('\n');

const linesTokens = (function () {
    const tokens: number[] = [];
    let last = 0;
    for (const l of lines) {
        tokens.push(last, last + l.length);
        last += l.length + 1;
    }
    if (tokens[tokens.length - 1] > linesData.length) tokens[tokens.length - 1] = linesData.length;
    return tokens;
}());

describe('fixed text column', () => {
    const col = FixedColumn({ data: linesData, indices: linesTokens, count: lines.length });
    const col1 = col(0, 5, Column.Schema.float);
    const col2 = col(5, 4, Column.Schema.str);
    it('number', () => {
        expect(col1.value(0)).toBe(1.123);
        expect(col1.value(1)).toBe(1.0);
        expect(col1.value(2)).toBe(1.1);
        expect(col1.value(3)).toBe(0);
        expect(col1.value(4)).toBe(5);
    });
    it('str', () => {
        expect(col2.value(0)).toBe('abc');
        expect(col2.value(1)).toBe('a');
        expect(col2.value(2)).toBe('bc');
        expect(col2.value(3)).toBe('');
        expect(col2.value(4)).toBe('');
    });
});

describe('token text column', () => {
    const tokensData = '321';
    const col = TokenColumn({ data: tokensData, indices: [0, 1, 1, 2, 2, 3], count: 3 });
    const col1 = col(Column.Schema.int);
    it('number', () => {
        expect(col1.value(0)).toBe(3);
        expect(col1.value(1)).toBe(2);
        expect(col1.value(2)).toBe(1);
    });
});

describe('binary column', () => {
    it('window works', () => {
        const xs = new Float64Array([1, 2, 3, 4]);
        const w1 = ColumnHelpers.typedArrayWindow(xs, { start: 1 });
        const w2 = ColumnHelpers.typedArrayWindow(xs, { start: 2, end: 4 });

        expect(w1.length).toBe(3);
        for (let i = 0; i < w1.length; i++) expect(w1[i]).toBe(xs[i + 1]);

        expect(w2.length).toBe(2);
        for (let i = 0; i < w2.length; i++) expect(w2[i]).toBe(xs[i + 2]);
    });
});
