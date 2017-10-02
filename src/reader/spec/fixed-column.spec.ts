/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import FixedColumn from '../common/text/column/fixed'
import { ColumnType } from '../common/column'

const lines = [
    '1.123 abc',
    '1.00  a',
    '1.1    bcd   ',
    '',
    ' 5'
]

const data = lines.join('\n');

const linesTokens = (function () {
    const tokens: number[] = [];
    let last = 0;
    for (const l of lines) {
        tokens.push(last, last + l.length);
        last += l.length + 1;
    }
    if (tokens[tokens.length - 1] > data.length) tokens[tokens.length - 1] = data.length;
    return tokens;
}());

describe('fixed text column', () => {
    const col = FixedColumn({ data, tokens: linesTokens, count: lines.length });
    const col1 = col(0, 5, ColumnType.float);
    const col2 = col(5, 4, ColumnType.str);
    it('number', () => {
        expect(col1.value(0)).toBe(1.123);
        expect(col1.value(1)).toBe(1.0);
        expect(col1.value(2)).toBe(1.1);
        expect(col1.value(3)).toBe(0);
        expect(col1.value(4)).toBe(5);
    })
    it('str', () => {
        expect(col2.value(0)).toBe('abc');
        expect(col2.value(1)).toBe('a');
        expect(col2.value(2)).toBe('bc');
        expect(col2.value(3)).toBe('');
        expect(col2.value(4)).toBe('');
    })
});
