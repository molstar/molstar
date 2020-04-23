/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as CCP4 from '../ccp4/parser';

function createCcp4Data() {
    const data = new Uint8Array(4 * 256 + 6);

    const dv = new DataView(data.buffer);

    dv.setInt8(52 * 4, 'M'.charCodeAt(0));
    dv.setInt8(52 * 4 + 1, 'A'.charCodeAt(0));
    dv.setInt8(52 * 4 + 2, 'P'.charCodeAt(0));
    dv.setInt8(52 * 4 + 3, ' '.charCodeAt(0));

    dv.setUint8(53 * 4, 17);
    dv.setUint8(53 * 4 + 1, 17);

    dv.setInt32(0 * 4, 1); // NC
    dv.setInt32(1 * 4, 2); // NR
    dv.setInt32(2 * 4, 3); // NS

    dv.setInt32(3 * 4, 0); // MODE

    return data;
}

describe('ccp4 reader', () => {
    it('basic', async () => {
        const data = createCcp4Data();
        const parsed = await CCP4.parse(data, 'test.ccp4').run();

        if (parsed.isError) {
            throw new Error(parsed.message);
        }

        const ccp4File = parsed.result;
        const { header } = ccp4File;

        expect(header.NC).toBe(1);
        expect(header.NR).toBe(2);
        expect(header.NS).toBe(3);
    });
});
