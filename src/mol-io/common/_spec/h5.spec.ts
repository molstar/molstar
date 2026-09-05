/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Smoke test for the vendored jsfive HDF5 reader.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { readFileSync } from 'fs';
import * as path from 'path';
import { File as H5File } from '../h5';

function loadFixture(name: string): ArrayBuffer {
    const buf = readFileSync(path.resolve(__dirname, 'fixtures', name));
    return buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength);
}

describe('h5 (jsfive) smoke', () => {
    let f: any;

    beforeAll(() => {
        const ab = loadFixture('h5-smoke.h5');
        f = new H5File(ab as any, 'h5-smoke.h5');
    });

    it('reads file-level attributes', () => {
        expect(f.attrs['note']).toBe('jsfive smoke test');
    });

    it('reads a scalar dataset in a nested group', () => {
        const pi = f.get('scalars/pi');
        expect(pi).toBeTruthy();
        // jsfive returns scalars wrapped in a length-1 array
        expect(pi.value[0]).toBeCloseTo(3.14159, 5);
    });

    it('reads a contiguous int dataset', () => {
        const v = f.get('vec_int').value;
        expect(Array.from(v as Iterable<number>)).toEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    });

    it('reads a gzip-compressed chunked float dataset', () => {
        const v = f.get('vec_gzip').value;
        expect(v.length).toBe(2048);
        expect(v[0]).toBe(0);
        expect(v[1024]).toBe(1024);
        expect(v[2047]).toBe(2047);
    });

    it('reads a 2D dataset and its attributes', () => {
        const m = f.get('a/b/m');
        expect(m.shape).toEqual([2, 3]);
        expect(Array.from(m.value as Iterable<number>)).toEqual([1, 2, 3, 4, 5, 6]);
        expect(m.attrs['units']).toBe('pixels');
    });
});
