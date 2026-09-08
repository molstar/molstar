/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Color } from '../../../mol-util/color';
import { bodyMaskFileName, maskBaseName, orderBodiesBySize, resolveBodyMaskParams } from '../internal/export';
import { BodyInfo } from '../types';

function body(id: number, voxelCount: number, extra: Partial<BodyInfo> = {}): BodyInfo {
    return { id, name: `Body ${id}`, color: Color(0xff0000), views: [], remainder: false, voxelCount, ...extra };
}

describe('export helpers', () => {
    it('orders bodies largest first and drops empty ones', () => {
        const ordered = orderBodiesBySize([body(1, 10), body(2, 50), body(3, 0), body(4, 50)]);
        expect(ordered.map(b => b.id)).toEqual([2, 4, 1]);
    });

    it('builds file names', () => {
        expect(maskBaseName('emd_1234.map.gz')).toBe('emd_1234');
        expect(maskBaseName('run_class001.mrc')).toBe('run_class001');
        expect(maskBaseName(undefined)).toBe('volume');
        expect(bodyMaskFileName('emd_1234', 1)).toBe('emd_1234_body001_mask.mrc');
        expect(bodyMaskFileName('x', 12)).toBe('x_body012_mask.mrc');
    });

    it('resolves per-body overrides', () => {
        const defaults = { extend: 6, softEdge: 2, pruneBelowThreshold: true };
        expect(resolveBodyMaskParams(body(1, 1), defaults)).toEqual(defaults);
        expect(resolveBodyMaskParams(body(1, 1, { extend: 1 }), defaults)).toEqual({ extend: 1, softEdge: 2, pruneBelowThreshold: true });
    });
});
