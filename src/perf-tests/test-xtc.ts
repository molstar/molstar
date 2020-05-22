/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';
import { parseXtc } from '../mol-io/reader/xtc/parser';

console.log('reading');
console.time('read');
fs.readFile('C:\\Projects\\mol-star\\molstar\\build\\tests\\test.xtc', async (err, data) => {
    console.log(err);
    console.timeEnd('read');
    console.time('parse');
    const ret = await parseXtc(new Uint8Array(data)).run(o => {
        console.log(`${o.root.progress.current}/${o.root.progress.max}`);
    }, 1000);
    console.timeEnd('parse');

    if (ret.isError) {
        console.log(ret.message);
    } else {
        console.log(ret.result?.frames.length);
        console.log(ret.result?.frames[0].x[250]);
    }
});