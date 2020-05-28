#!/usr/bin/env node
/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse';
import * as fs from 'fs';
import * as path from 'path';
import fetch from 'node-fetch';
import { UniqueArray } from '../../mol-data/generic';

const LIPIDS_DIR = path.resolve(__dirname, '../../../../build/lipids/');

const MARTINI_LIPIDS_PATH = path.resolve(LIPIDS_DIR, 'martini_lipids.itp');
const MARTINI_LIPIDS_URL = 'http://www.cgmartini.nl/images/parameters/lipids/Collections/martini_v2.0_lipids_all_201506.itp';

async function ensureAvailable(path: string, url: string) {
    if (FORCE_DOWNLOAD || !fs.existsSync(path)) {
        const name = url.substr(url.lastIndexOf('/') + 1);
        console.log(`downloading ${name}...`);
        const data = await fetch(url);
        if (!fs.existsSync(LIPIDS_DIR)) {
            fs.mkdirSync(LIPIDS_DIR);
        }
        fs.writeFileSync(path, await data.text());
        console.log(`done downloading ${name}`);
    }
}

async function ensureLipidsAvailable() { await ensureAvailable(MARTINI_LIPIDS_PATH, MARTINI_LIPIDS_URL); }

async function run(out: string) {
    await ensureLipidsAvailable();
    const lipidsItpStr = fs.readFileSync(MARTINI_LIPIDS_PATH, 'utf8');

    const lipids = UniqueArray.create<string>();
    const reLipid = /\[moleculetype\]\n; molname      nrexcl\n +([a-zA-Z]{3,5})/g;
    let m: RegExpExecArray | null;

    while ((m = reLipid.exec(lipidsItpStr)) !== null) {
        const v = m[0].substr(m[0].lastIndexOf(' ') + 1);
        UniqueArray.add(lipids, v, v);
    }


    const lipidNames = JSON.stringify(lipids.array);

    if (out) {
        const output = `/**
* Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
*
* Code-generated lipid params file. Names extracted from Martini FF lipids itp.
*
* @author molstar/lipid-params cli
*/

export const LipidNames = new Set(${lipidNames.replace(/"/g, "'").replace(/,/g, ', ')});
`;
        fs.writeFileSync(out, output);
    } else {
        console.log(lipidNames);
    }
}

const parser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Create lipid params (from martini lipids itp)'
});
parser.addArgument([ '--out', '-o' ], {
    help: 'Generated lipid params output path, if not given printed to stdout'
});
parser.addArgument([ '--forceDownload', '-f' ], {
    action: 'storeTrue',
    help: 'Force download of martini lipids itp'
});
interface Args {
    out: string
    forceDownload: boolean
}
const args: Args = parser.parseArgs();

const FORCE_DOWNLOAD = args.forceDownload;

run(args.out || '').catch(e => {
    console.error(e);
});