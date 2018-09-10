/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as argparse from 'argparse'
import { preprocessFile } from './preprocess/preprocess';

const cmdParser = new argparse.ArgumentParser({
    addHelp: true,
    description: 'Preprocess CIF files to include custom properties and convert them to BinaryCIF format.'
});
cmdParser.addArgument(['--input', '-i'], { help: 'Input filename', required: true });
cmdParser.addArgument(['--outCIF', '-oc'], { help: 'Output CIF filename', required: false });
cmdParser.addArgument(['--outBCIF', '-ob'], { help: 'Output BinaryCIF filename', required: false });

// TODO: "bulk" mode

interface CmdArgs {
    input: string,
    outCIF?: string,
    outBCIF?: string
}

const cmdArgs = cmdParser.parseArgs() as CmdArgs;

if (cmdArgs.input) preprocessFile(cmdArgs.input, cmdArgs.outCIF, cmdArgs.outBCIF);

// example:
// node build\node_modules\servers\model\preprocess -i e:\test\Quick\1cbs_updated.cif -oc e:\test\mol-star\model\1cbs.cif -ob e:\test\mol-star\model\1cbs.bcif