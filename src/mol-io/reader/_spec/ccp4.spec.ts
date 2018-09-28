/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import CCP4 from '../ccp4/parser'
import { FileHandle } from '../../common/file-handle';

const ccp4Buffer = new Uint8Array([])

describe('ccp4 reader', () => {
    it('basic', async () => {
        const file = FileHandle.fromBuffer(ccp4Buffer)
        const parsed = await CCP4(file).run();

        if (parsed.isError) {
            console.log(parsed)
            return;
        }
        // const ccp4File = parsed.result;
        // const { header, values } = ccp4File;

        // TODO
    });
});
