#!/usr/bin/env node
/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse';
import * as util from 'util';
import * as fs from 'fs';
import * as zlib from 'zlib';
import convert from './converter';

require('util.promisify').shim();

async function process(srcPath: string, outPath: string, configPath?: string, filterPath?: string) {
    const config = configPath ? JSON.parse(fs.readFileSync(configPath, 'utf8')) : void 0;
    const filter = filterPath ? fs.readFileSync(filterPath, 'utf8') : void 0;

    const res = await convert(srcPath, false, config, filter);
    await write(outPath, res);
}

const zipAsync = util.promisify<zlib.InputType, Buffer>(zlib.gzip);

async function write(outPath: string, res: Uint8Array) {
    const isGz = /\.gz$/i.test(outPath);
    if (isGz) {
        res = await zipAsync(res);
    }
    fs.writeFileSync(outPath, res);
}

function run(args: Args) {
    process(args.src, args.out, args.config, args.filter);
}

const parser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Convert any CIF file to a BCIF file'
});
parser.addArgument([ 'src' ], {
    help: 'Source CIF path'
});
parser.addArgument([ 'out' ], {
    help: 'Output BCIF path'
});
parser.addArgument([ '-c', '--config' ], {
    help: 'Optional encoding strategy/precision config path',
    required: false
});
parser.addArgument([ '-f', '--filter' ], {
    help: 'Optional filter whitelist/blacklist path',
    required: false
});

interface Args {
    src: string
    out: string
    config?: string
    filter?: string
}
const args: Args = parser.parseArgs();

if (args) {
    run(args);
}