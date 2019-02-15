/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as CCP4 from '../ccp4/parser'

const ccp4Buffer = new Uint8Array(4 * 64)

describe('ccp4 reader', () => {
    it('basic', async () => {
        const parsed = await CCP4.parse(ccp4Buffer).run();

        if (parsed.isError) {
            console.log(parsed)
            return;
        }
        // const ccp4File = parsed.result;
        // const { header, values } = ccp4File;

        // TODO
    });
});
