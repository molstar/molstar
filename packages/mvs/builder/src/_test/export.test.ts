import { unzipSync } from 'fflate';
import { createMVSX } from '../export.js';
import { MVSData } from '../mvs-data.js';

test('MVSX output is deterministic across different clock times with explicit mtime', async () => {
    jest.useFakeTimers();
    try {
        const data = MVSData.createBuilder().getState();
        const assets = [{ name: 'asset.txt', content: 'archive asset' }];
        // ZIP stores local date fields; use matching fields for reproducibility across time zones.
        const options = { zip: { mtime: new Date(1980, 0, 1), level: 0 as const } };
        jest.setSystemTime(new Date('2025-01-01T12:00:00Z'));
        const first = await createMVSX(data, assets, options);
        jest.setSystemTime(new Date('2026-10-04T18:00:00Z'));
        const second = await createMVSX(data, assets, options);
        expect(second).toEqual(first);
        const files = unzipSync(first);
        expect(new TextDecoder().decode(files['asset.txt'])).toBe('archive asset');
        expect(JSON.parse(new TextDecoder().decode(files['index.mvsj']))).toEqual(data);
    } finally {
        jest.useRealTimers();
    }
});

test('MVSX output uses the current time when ZIP options are omitted', async () => {
    jest.useFakeTimers();
    try {
        const data = MVSData.createBuilder().getState();
        jest.setSystemTime(new Date('2025-01-01T12:00:00Z'));
        const first = await createMVSX(data, []);
        jest.setSystemTime(new Date('2026-10-04T18:00:00Z'));
        const second = await createMVSX(data, []);
        expect(second).not.toEqual(first);
        expect(unzipSync(second)).toEqual(unzipSync(first));
    } finally {
        jest.useRealTimers();
    }
});
