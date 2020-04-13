/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';
import * as path from 'path';
import * as argparse from 'argparse';
import { runMaster, PreprocessEntry } from './parallel';
import { ModelPropertyProviderConfig } from '../property-provider';

function description() {
    const exampleCfg = {
        numProcesses: 1,
        customProperties: {
            sources: [
                'wwpdb'
            ],
            params: {
                wwPDB: {
                    chemCompBondTablePath: './build/data/ccb.bcif'
                }
            }
        }
    };

    return `Preprocess CIF files to include custom properties and convert them to BinaryCIF format.\n\nExample cfg.json: ${JSON.stringify(exampleCfg, null, 2)}`;
}

const cmdParser = new argparse.ArgumentParser({
    addHelp: true,
    description: description()
});
cmdParser.addArgument(['--input', '-i'], { help: 'Input filename', required: false });
cmdParser.addArgument(['--outCIF', '-oc'], { help: 'Output CIF filename', required: false });
cmdParser.addArgument(['--outBCIF', '-ob'], { help: 'Output BinaryCIF filename', required: false });
// TODO: add back? cmdParser.addArgument(['--bulk', '-b'], { help: 'Bulk JSON ({ numProcesses?: number, entries: { source: string, cif?: string, bcif?: string }[] })', required: false });
cmdParser.addArgument(['--cfg', '-c'], { help: 'Config file path', required: false });
cmdParser.addArgument(['--folderIn', '-fin'], { help: 'Convert folder', required: false });
cmdParser.addArgument(['--folderOutCIF', '-foc'], { help: 'Convert folder text output', required: false });
cmdParser.addArgument(['--folderOutBCIF', '-fob'], { help: 'Convert folder binary output', required: false });
cmdParser.addArgument(['--folderNumProcesses', '-fp'], { help: 'Convert folder num processes', required: false });

interface CmdArgs {
    // bulk?: string,
    help?: any,
    cfg?: string,
    input?: string,
    outCIF?: string,
    outBCIF?: string,
    folderIn?: string,
    folderOutCIF?: string,
    folderOutBCIF?: string,
    folderNumProcesses?: string
}


export interface PreprocessConfig {
    numProcesses?: number,
    customProperties?: ModelPropertyProviderConfig | string
}

const cmdArgs = cmdParser.parseArgs() as CmdArgs;

if (Object.keys(cmdArgs).filter(k => (cmdArgs as any)[k] !== null).length === 0 || typeof cmdArgs.help !== 'undefined') {
    cmdParser.printHelp();
    process.exit(0);
}

let entries: PreprocessEntry[] = [];
let config: PreprocessConfig = { numProcesses: cmdArgs.folderIn ? +(cmdArgs.folderNumProcesses || 1) : 1, customProperties: void 0 };

if (cmdArgs.input) entries.push({ source: cmdArgs.input, cif: cmdArgs.outCIF, bcif: cmdArgs.outBCIF });
// else if (cmdArgs.bulk) runBulk(cmdArgs.bulk);
else if (cmdArgs.folderIn) findEntries();

if (cmdArgs.cfg) {
    config = JSON.parse(fs.readFileSync(cmdArgs.cfg, 'utf8')) as PreprocessConfig;
}

runMaster(config, entries);

function findEntries() {
    const files = fs.readdirSync(cmdArgs.folderIn!);
    const cifTest = /\.cif$/;
    for (const f of files) {
        if (!cifTest.test(f)) continue;

        entries.push({
            source: path.join(cmdArgs.folderIn!, f),
            cif: cmdArgs.folderOutCIF ? path.join(cmdArgs.folderOutCIF!, f) : void 0,
            bcif: cmdArgs.folderOutBCIF ? path.join(cmdArgs.folderOutBCIF!, path.parse(f).name + '.bcif') : void 0,
        });
    }
}

// example:
// node build\node_modules\servers\model\preprocess -i e:\test\Quick\1cbs_updated.cif -oc e:\test\mol-star\model\1cbs.cif -ob e:\test\mol-star\model\1cbs.bcif