#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 *
 * Command-line application for converting MolViewSpec MVSJ into MSVX files
 * Build: npm run build
 * Run:   node lib/commonjs/cli/mvs/mvs-mvsj-to-mvsx -i examples/mvs/1cbs.mvsj -o tmp/1cbs.mvsx
 */

import { ArgumentParser } from 'argparse';
import fs from 'fs';
import { MVSData } from '../../extensions/mvs/mvs-data';
import { setFSModule } from '../../mol-util/data-source';


setFSModule(fs);

/** Command line argument values for `main` */
interface Args {
    input: string[],
    output: string[] | undefined,
    base_uri: string | undefined,
    skip_external: boolean,
}

/** Return parsed command line arguments for `main` */
function parseArguments(): Args {
    const parser = new ArgumentParser({ description: 'Command-line application for converting MolViewSpec MVSJ into MSVX files' });
    parser.add_argument('-i', '--input', { required: true, nargs: '+', help: 'Input file(s) in .mvsj format.' });
    parser.add_argument('-o', '--output', { required: false, nargs: '+', help: 'File path(s) for output files in .mvsx format (one output path for each input file). If ommitted, filenames will be created automatically by replacing file extension.' });
    parser.add_argument('--base-uri', { help: 'Base URI/path used to resolve relative URIs in the input file (default: path of the input file itself). Use `--base-uri .` for using the current working directory as base URI.' });
    parser.add_argument('--skip-external', { action: 'store_true', help: 'Do not include external resources (i.e. absolute URIs) in the MVSX.' });
    const args: Args = parser.parse_args();
    if (args.output && args.output.length !== args.input.length) {
        parser.error(`argument: --output: must specify the same number of input and output file paths (specified ${args.input.length} input path${args.input.length !== 1 ? 's' : ''} but ${args.output.length} output path${args.output.length !== 1 ? 's' : ''})`);
    }
    return { ...args };
}

/** Main workflow for converting MVSJ to MVSX files. */
async function main(args: Args): Promise<void> {
    const cache = {};
    for (let i = 0; i < args.input.length; i++) {
        const input = args.input[i];
        const output = args.output?.[i] ?? input.replace(/(\.mvsj)?$/i, '.mvsx');
        console.log(`Processing ${input} -> ${output}`);
        const mvsj = fs.readFileSync(input, { encoding: 'utf8' });
        const mvsData = MVSData.fromMVSJ(mvsj);
        const mvsx = await MVSData.toMVSX(mvsData, {
            baseUri: args.base_uri ?? input,
            skipExternal: args.skip_external,
            cache,
        });
        fs.writeFileSync(output, mvsx);
    }
}

main(parseArguments());
