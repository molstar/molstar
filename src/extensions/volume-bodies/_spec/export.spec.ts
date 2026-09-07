/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { parse } from '../../../mol-io/reader/ccp4/parser';
import { CCP4Writer } from '../../../mol-io/writer/ccp4/ccp4';
import { volumeFromCcp4 } from '../../../mol-model-formats/volume/ccp4';
import { Color } from '../../../mol-util/color';
import { bodyMaskFileName, maskBaseName, orderBodiesBySize, resolveBodyMaskParams } from '../internal/export';
import { BodyInfo } from '../types';
import { CanonicalOrder, SwappedOrder, createTestVolume } from './test-volume';

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

describe('CCP4Writer round trip', () => {
    for (const order of [CanonicalOrder, SwappedOrder]) {
        it(`preserves voxel values for axis order ${order}`, async () => {
            const volume = createTestVolume([4, 3, 2], (i, j, k) => i + 10 * j + 100 * k, order);
            const data = volume.grid.cells.data as Float32Array;
            const buffer = CCP4Writer.writeMrc(volume.grid, data);
            expect(buffer.byteLength).toBe(1024 + 4 * 3 * 2 * 4);

            const parsed = await parse(new Uint8Array(buffer), 'test.mrc').run();
            if (parsed.isError) throw new Error(parsed.message);
            expect(parsed.result.header.NC).toBe(4);
            expect(parsed.result.header.NR).toBe(3);
            expect(parsed.result.header.NS).toBe(2);
            expect(parsed.result.header.MODE).toBe(2);
            expect(parsed.result.header.AMIN).toBe(0);
            expect(parsed.result.header.AMAX).toBe(123);

            const read = await volumeFromCcp4(parsed.result).run();
            const { space, data: values } = read.grid.cells;
            expect(Array.from(space.dimensions)).toEqual([4, 3, 2]);
            for (let i = 0; i < 4; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 2; k++) {
                expect(space.get(values, i, j, k)).toBe(i + 10 * j + 100 * k);
            }
        });
    }
});
