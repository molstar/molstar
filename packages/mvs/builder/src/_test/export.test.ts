import { unzipSync } from 'fflate';
import { createMVSX } from '../export.js';
import { MVSData } from '../mvs-data.js';

test('MVSX output is deterministic across different clock times', async () => {
    jest.useFakeTimers();
    try {
        const data = MVSData.createBuilder().getState();
        const assets = [{ name: 'asset.txt', content: 'archive asset' }];
        jest.setSystemTime(new Date('2025-01-01T12:00:00Z'));
        const first = await createMVSX(data, assets);
        jest.setSystemTime(new Date('2026-10-04T18:00:00Z'));
        const second = await createMVSX(data, assets);
        expect(second).toEqual(first);
        const files = unzipSync(first);
        expect(new TextDecoder().decode(files['asset.txt'])).toBe('archive asset');
        expect(JSON.parse(new TextDecoder().decode(files['index.mvsj']))).toEqual(data);
    } finally {
        jest.useRealTimers();
    }
});
