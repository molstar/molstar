/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { AnnotationRow } from '../schemas';
import { groupRows } from '../selections';
import { range } from '../utils';


describe('groupRows', () => {
    it('groupRows', async () => {
        const rows = [
            { label: 'A' }, { label: 'B', group_id: 1 }, { label: 'C', group_id: 'x' }, { label: 'D', group_id: 1 },
            { label: 'E' }, { label: 'F' }, { label: 'G', group_id: 'x' }, { label: 'H', group_id: 'x' },
        ] as any as AnnotationRow[];
        const g = groupRows(rows);
        const groupedIndices = range(g.count).map(i => g.grouped.slice(g.offsets[i], g.offsets[i + 1]));
        const groupedRows = groupedIndices.map(group => group.map(j => rows[j]));
        expect(groupedRows).toEqual([
            [{ label: 'A' }],
            [{ label: 'B', group_id: 1 }, { label: 'D', group_id: 1 }],
            [{ label: 'C', group_id: 'x' }, { label: 'G', group_id: 'x' }, { label: 'H', group_id: 'x' }],
            [{ label: 'E' }],
            [{ label: 'F' }],
        ]);
    });
});
