#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse';
import * as fs from 'fs';
import * as path from 'path';

import { parseSyminfoLib } from './util/parse';
import { generateSpacegroupSyminfoTs, generateSyminfoLibSpecTs } from './util/generate';

const DATA_PATH = path.resolve(__dirname, '../../../../data/sym/syminfo.lib');
const SYMINFO_TS_PATH = path.resolve(__dirname, '../../../../src/mol-math/geometry/spacegroup/syminfo.ts');
const SYMINFO_LIB_SPEC_TS_PATH = path.resolve(__dirname, '../../../../src/mol-math/geometry/spacegroup/_spec/syminfo.lib.ts');

const parser = new argparse.ArgumentParser({
    add_help: true,
    description: 'Generate mol* spacegroup data from CCP4 syminfo.lib (data/sym/syminfo.lib)'
});
parser.add_argument('--dataPath', '-d', { default: DATA_PATH, help: 'Path to the raw syminfo.lib text' });
parser.add_argument('--out', '-o', { default: SYMINFO_TS_PATH, help: 'Output path for src/mol-math/geometry/spacegroup/syminfo.ts (RawSpacegroupData)' });
parser.add_argument('--outSpec', '-os', { default: SYMINFO_LIB_SPEC_TS_PATH, help: 'Output path for _spec/syminfo.lib.ts (the SyminfoEntry[] regression oracle)' });

interface Args {
    dataPath: string
    out: string
    outSpec: string
}

function run() {
    const args: Args = parser.parse_args();

    const raw = fs.readFileSync(args.dataPath, 'utf8');
    const entries = parseSyminfoLib(raw);
    console.log(`parsed ${entries.length} settings from '${args.dataPath}'`);

    fs.writeFileSync(args.out, generateSpacegroupSyminfoTs(entries));
    console.log(`wrote '${args.out}'`);

    fs.writeFileSync(args.outSpec, generateSyminfoLibSpecTs(entries));
    console.log(`wrote '${args.outSpec}'`);
}

run();
