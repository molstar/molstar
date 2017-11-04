/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs'
import convert from './converter'

(async function () {
    if (process.argv.length !== 4) {
        console.log('Usage:\nnode cif2bcif input.cif output.bcif');
        return;
    }
    const src = process.argv[2];
    const out = process.argv[3];

    const res = await convert(src);
    fs.writeFileSync(out, res);
}());