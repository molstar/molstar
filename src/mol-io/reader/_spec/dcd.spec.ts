/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseDcd } from '../dcd/parser';

function createDcdData() {
    const data = new Uint8Array(4 * 128);

    const dv = new DataView(data.buffer);

    // set little endian
    dv.setInt32(0, 84);

    // set format string
    dv.setUint8(4, 'D'.charCodeAt(0));
    dv.setUint8(5, 'R'.charCodeAt(0));
    dv.setUint8(6, 'O'.charCodeAt(0));
    dv.setUint8(7, 'C'.charCodeAt(0));

    dv.setInt32(8, 1); // NSET

    // header end
    dv.setInt32(22 * 4, 84);

    // title
    const titleEnd = 164;
    const titleStart = 23 * 4 + 1;
    dv.setInt32(23 * 4, titleEnd);
    dv.setInt32(titleStart + titleEnd + 4 - 1, titleEnd);

    // atoms
    const atomStart = 23 * 4 + titleEnd + 8;
    dv.setInt32(atomStart, 4);
    dv.setInt32(atomStart + 4, 1); // one atom
    dv.setInt32(atomStart + 8, 4);

    // coords
    const coordsStart = atomStart + 12;
    dv.setInt32(coordsStart, 4);
    dv.setFloat32(coordsStart + 4, 0.1);
    dv.setInt32(coordsStart + 8, 4);
    dv.setInt32(coordsStart + 12, 4);
    dv.setFloat32(coordsStart + 16, 0.2);
    dv.setInt32(coordsStart + 20, 4);
    dv.setInt32(coordsStart + 24, 4);
    dv.setFloat32(coordsStart + 28, 0.3);
    dv.setInt32(coordsStart + 32, 4);

    return data;
}

describe('dcd reader', () => {
    it('basic', async () => {
        const data = createDcdData();
        const parsed = await parseDcd(data).run();

        if (parsed.isError) {
            throw new Error(parsed.message);
        }

        const dcdFile = parsed.result;
        const { header, frames } = dcdFile;

        expect(header.NSET).toBe(1);
        expect(header.NATOM).toBe(1);

        expect(frames[0].x[0]).toBeCloseTo(0.1, 0.0001);
        expect(frames[0].y[0]).toBeCloseTo(0.2, 0.0001);
        expect(frames[0].z[0]).toBeCloseTo(0.3, 0.0001);
    });
});
