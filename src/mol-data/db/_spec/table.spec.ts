/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as ColumnHelpers from '../column-helpers';
import Column from '../column';
import Table from '../table';

describe('column', () => {
    const cc = Column.ofConst(10, 2, Column.Schema.int);
    const arr = Column.ofArray({ array: [1, 2, 3, 4], schema: Column.Schema.int });
    const arrNumberList = Column.ofArray({ array: [[1, 2], [3, 4], [5, 6]], schema: Column.Schema.List(' ', x => parseInt(x, 10)) });
    const arrStringList = Column.ofArray({ array: [['a', 'b'], ['c', 'd'], ['e', 'f']], schema: Column.Schema.List(',', x => x) });
    const arrWindow = Column.window(arr, 1, 3);

    const typed = Column.ofArray({ array: new Int32Array([1, 2, 3, 4]), schema: Column.Schema.int });
    const typedWindow = Column.window(typed, 1, 3);

    const numStr = Column.ofArray({ array: [1, 2] as any, schema: Column.Schema.str });

    it('constant', () => {
        expect(cc.rowCount).toBe(2);
        expect(cc.value(0)).toBe(10);
    });

    it('arr', () => {
        expect(arr.rowCount).toBe(4);
        expect(arr.value(1)).toBe(2);
        expect(arrWindow.value(0)).toBe(2);
        expect(arrWindow.rowCount).toBe(2);
    });

    it('arrList', () => {
        expect(arrNumberList.rowCount).toBe(3);
        expect(arrNumberList.value(1)).toEqual([3, 4]);
        expect(arrStringList.rowCount).toBe(3);
        expect(arrStringList.value(2)).toEqual(['e', 'f']);
    });

    it('typed', () => {
        expect(typedWindow.value(0)).toBe(2);
        expect(typedWindow.rowCount).toBe(2);
        expect(ColumnHelpers.isTypedArray(typedWindow.toArray())).toBe(true);
    });

    it('numStr', () => {
        expect(numStr.value(0)).toBe('1');
        expect(numStr.toArray()).toEqual(['1', '2']);
    });

    it('view', () => {
        expect(Column.view(arr, [1, 0, 3, 2]).toArray()).toEqual([2, 1, 4, 3]);
        expect(Column.view(arr, [1, 3]).toArray()).toEqual([2, 4]);
    });

    it('map to array', () => {
        expect(Column.mapToArray(arrWindow, x => x + 1)).toEqual([3, 4]);
    });
});

describe('table', () => {
    const schema = {
        x: Column.Schema.int,
        n: Column.Schema.str
    };

    it('ofRows', () => {
        const t = Table.ofRows(schema, [
            { x: 10, n: 'row1' },
            { x: -1, n: 'row2' },
        ]);
        expect(t.x.toArray()).toEqual([10, -1]);
        expect(t.n.toArray()).toEqual(['row1', 'row2']);
    });

    it('ofColumns', () => {
        const t = Table.ofColumns(schema, {
            x: Column.ofArray({ array: [10, -1], schema: Column.Schema.int }),
            n: Column.ofArray({ array: ['row1', 'row2'], schema: Column.Schema.str }),
        });
        expect(t.x.toArray()).toEqual([10, -1]);
        expect(t.n.toArray()).toEqual(['row1', 'row2']);
    });

    it('ofArrays', () => {
        const t = Table.ofArrays(schema, {
            x: [10, -1],
            n: ['row1', 'row2'],
        });
        expect(t.x.toArray()).toEqual([10, -1]);
        expect(t.n.toArray()).toEqual(['row1', 'row2']);
    });

    it('pickColumns', () => {
        const t = Table.ofColumns(schema, {
            x: Column.ofArray({ array: [10, -1], schema: Column.Schema.int }),
            n: Column.ofArray({ array: ['row1', 'row2'], schema: Column.Schema.str }),
        });
        const s = { x: Column.Schema.int, y: Column.Schema.int };
        const picked = Table.pickColumns(s, t, { y: Column.ofArray({ array: [3, 4], schema: Column.Schema.int })});
        expect(picked._columns).toEqual(['x', 'y']);
        expect(picked._rowCount).toEqual(2);
        expect(picked.x.toArray()).toEqual([10, -1]);
        expect(picked.y.toArray()).toEqual([3, 4]);
    });

    it('view', () => {
        const t = Table.ofColumns(schema, {
            x: Column.ofArray({ array: [10, -1], schema: Column.Schema.int }),
            n: Column.ofArray({ array: ['row1', 'row2'], schema: Column.Schema.str }),
        });
        const s = { x: Column.Schema.int };
        const view = Table.view(t, s, [1]);
        expect(view._columns).toEqual(['x']);
        expect(view._rowCount).toEqual(1);
        expect(view.x.toArray()).toEqual([-1]);
    });

    it('sort', () => {
        const t = Table.ofColumns<typeof schema>(schema, {
            x: Column.ofArray({ array: [10, -1], schema: Column.Schema.int }),
            n: Column.ofArray({ array: ['row1', 'row2'], schema: Column.Schema.str }),
        });
        const { x } = t;
        const sorted = Table.sort(t, (i, j) => x.value(i) - x.value(j));
        expect(sorted.x.toArray()).toEqual([-1, 10]);
        expect(sorted.n.toArray()).toEqual(['row2', 'row1']);
    });
});