/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { parseXtc } from '../../mol-io/reader/xtc/parser';
import './index.html';

const parent = document.getElementById('app')!;
const btn = document.createElement('button');
btn.innerText = 'run';
btn.onclick = run;
parent.appendChild(btn);

async function run() {
    const req = await fetch('test.xtc');
    const data = await req.arrayBuffer();
    console.log(data.byteLength);
    console.time('parse');
    const ret = await parseXtc(new Uint8Array(data)).run(o => {
        console.log(o.root.progress.current, o.root.progress.max);
    }, 1000);
    console.timeEnd('parse');
    console.log(ret);
    btn.innerText = 'done';
}