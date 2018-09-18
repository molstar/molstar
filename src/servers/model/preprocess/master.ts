/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs'
import * as path from 'path'
import * as argparse from 'argparse'
import { preprocessFile } from './preprocess';
import { ParallelPreprocessConfig, runMaster } from './parallel';

const cmdParser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Preprocess CIF files to include custom properties and convert them to BinaryCIF format.'
});
cmdParser.addArgument(['--input', '-i'], { help: 'Input filename', required: false });
cmdParser.addArgument(['--outCIF', '-oc'], { help: 'Output CIF filename', required: false });
cmdParser.addArgument(['--outBCIF', '-ob'], { help: 'Output BinaryCIF filename', required: false });
cmdParser.addArgument(['--bulk', '-b'], { help: 'Bulk JSON ({ numProcesses?: number, entries: { source: string, cif?: string, bcif?: string }[] })', required: false });
cmdParser.addArgument(['--folderIn', '-f'], { help: 'Convert folder', required: false });
cmdParser.addArgument(['--folderOutCIF', '-foc'], { help: 'Convert folder text output', required: false });
cmdParser.addArgument(['--folderOutBCIF', '-fob'], { help: 'Convert folder binary output', required: false });
cmdParser.addArgument(['--folderNumProcesses', '-fp'], { help: 'Convert folder num processes', required: false });

interface CmdArgs {
    bulk?: string,
    input?: string,
    outCIF?: string,
    outBCIF?: string,
    folderIn?: string,
    folderOutCIF?: string,
    folderOutBCIF?: string,
    folderNumProcesses?: string
}

const cmdArgs = cmdParser.parseArgs() as CmdArgs;

if (cmdArgs.input) preprocessFile(cmdArgs.input, cmdArgs.outCIF, cmdArgs.outBCIF);
else if (cmdArgs.bulk) runBulk(cmdArgs.bulk);
else if (cmdArgs.folderIn) runFolder(cmdArgs);

function runBulk(input: string) {
    const config = JSON.parse(fs.readFileSync(input, 'utf8')) as ParallelPreprocessConfig;
    runMaster(config);
}

function runFolder(args: CmdArgs) {
    const files = fs.readdirSync(args.folderIn!);
    const config: ParallelPreprocessConfig = { numProcesses: +args.folderNumProcesses! || 1, entries: [] };
    const cifTest = /\.cif$/;
    for (const f of files) {
        if (!cifTest.test(f)) continue;

        config.entries.push({
            source: path.join(args.folderIn!, f),
            cif: cmdArgs.folderOutCIF ? path.join(args.folderOutCIF!, f) : void 0,
            bcif: cmdArgs.folderOutBCIF ? path.join(args.folderOutBCIF!, path.parse(f).name + '.bcif') : void 0,
        });
    }
    runMaster(config);
}

// example:
// node build\node_modules\servers\model\preprocess -i e:\test\Quick\1cbs_updated.cif -oc e:\test\mol-star\model\1cbs.cif -ob e:\test\mol-star\model\1cbs.bcif