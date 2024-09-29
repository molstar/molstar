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

const MARTINI_LIPIDS_PATH = path.resolve(LIPIDS_DIR, 'martini_lipids_v3.itp');
const MARTINI_LIPIDS_URL = 'https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v3.0.0_phospholipids_v1.itp';

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

const extraLipids = ['DMPC'];
const v2lipids = ['DAPC', 'DBPC', 'DFPC', 'DGPC', 'DIPC', 'DLPC', 'DNPC', 'DOPC', 'DPPC', 'DRPC', 'DTPC', 'DVPC', 'DXPC', 'DYPC', 'LPPC', 'PAPC', 'PEPC', 'PGPC', 'PIPC', 'POPC', 'PRPC', 'PUPC', 'DAPE', 'DBPE', 'DFPE', 'DGPE', 'DIPE', 'DLPE', 'DNPE', 'DOPE', 'DPPE', 'DRPE', 'DTPE', 'DUPE', 'DVPE', 'DXPE', 'DYPE', 'LPPE', 'PAPE', 'PGPE', 'PIPE', 'POPE', 'PQPE', 'PRPE', 'PUPE', 'DAPS', 'DBPS', 'DFPS', 'DGPS', 'DIPS', 'DLPS', 'DNPS', 'DOPS', 'DPPS', 'DRPS', 'DTPS', 'DUPS', 'DVPS', 'DXPS', 'DYPS', 'LPPS', 'PAPS', 'PGPS', 'PIPS', 'POPS', 'PQPS', 'PRPS', 'PUPS', 'DAPG', 'DBPG', 'DFPG', 'DGPG', 'DIPG', 'DLPG', 'DNPG', 'DOPG', 'DPPG', 'DRPG', 'DTPG', 'DVPG', 'DXPG', 'DYPG', 'LPPG', 'PAPG', 'PGPG', 'PIPG', 'POPG', 'PRPG', 'DAPA', 'DBPA', 'DFPA', 'DGPA', 'DIPA', 'DLPA', 'DNPA', 'DOPA', 'DPPA', 'DRPA', 'DTPA', 'DVPA', 'DXPA', 'DYPA', 'LPPA', 'PAPA', 'PGPA', 'PIPA', 'POPA', 'PRPA', 'PUPA', 'DPP', 'DPPI', 'PAPI', 'PIPI', 'POP', 'POPI', 'PUPI', 'PVP', 'PVPI', 'PADG', 'PIDG', 'PODG', 'PUDG', 'PVDG', 'APC', 'CPC', 'IPC', 'LPC', 'OPC', 'PPC', 'TPC', 'UPC', 'VPC', 'BNSM', 'DBSM', 'DPSM', 'DXSM', 'PGSM', 'PNSM', 'POSM', 'PVSM', 'XNSM', 'DPCE', 'DXCE', 'PNCE', 'XNCE'];

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

    for (const v of extraLipids) {
        UniqueArray.add(lipids, v, v);
    }

    for (const v of v2lipids) {
        UniqueArray.add(lipids, v, v);
    }

    const lipidNames = JSON.stringify(lipids.array);

    if (out) {
        const output = `/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
    add_help: true,
    description: 'Create lipid params (from martini lipids itp)'
});
parser.add_argument('--out', '-o', {
    help: 'Generated lipid params output path, if not given printed to stdout'
});
parser.add_argument('--forceDownload', '-f', {
    action: 'store_true',
    help: 'Force download of martini lipids itp'
});
interface Args {
    out: string
    forceDownload: boolean
}
const args: Args = parser.parse_args();

const FORCE_DOWNLOAD = args.forceDownload;

run(args.out || '').catch(e => {
    console.error(e);
});