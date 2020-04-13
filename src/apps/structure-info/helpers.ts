/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as util from 'util';
import * as fs from 'fs';
import fetch from 'node-fetch';
require('util.promisify').shim();

import { CIF } from '../../mol-io/reader/cif';
import { Progress } from '../../mol-task';

const readFileAsync = util.promisify(fs.readFile);

async function readFile(path: string) {
    if (path.match(/\.bcif$/)) {
        const input = await readFileAsync(path);
        const data = new Uint8Array(input.byteLength);
        for (let i = 0; i < input.byteLength; i++) data[i] = input[i];
        return data;
    } else {
        return readFileAsync(path, 'utf8');
    }
}

async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data);
    const parsed = await comp.run(p => console.log(Progress.format(p)), 250);
    if (parsed.isError) throw parsed;
    return parsed.result;
}

export async function openCif(path: string) {
    const data = await readFile(path);
    return parseCif(data);
}

export async function downloadCif(url: string, isBinary: boolean) {
    const data = await fetch(url);
    return parseCif(isBinary ? new Uint8Array(await data.arrayBuffer()) : await data.text());
}