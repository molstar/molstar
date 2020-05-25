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

const BUILD_DIR = path.resolve(__dirname, '../build/');
const LIPIDS_DIR = path.resolve(BUILD_DIR, 'lipids/');

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

async function run() {
    await ensureLipidsAvailable();
    const lipidsItpStr = fs.readFileSync(MARTINI_LIPIDS_PATH, 'utf8');

    const m = lipidsItpStr.match(/\[moleculetype\]\n; molname      nrexcl\n(DGPC)/g);
    console.log(m);
}

const parser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Create lipid params (from martini lipids itp)'
});
parser.addArgument([ '--forceDownload', '-f' ], {
    action: 'storeTrue',
    help: 'Force download of martini lipids itp'
});
interface Args {
    forceDownload: boolean
}
const args: Args = parser.parseArgs();

const FORCE_DOWNLOAD = args.forceDownload;

run().catch(e => {
    console.error(e);
});