#!/usr/bin/env node
/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 *
 * Command-line application for validating MolViewSpec files
 * Build: npm run build
 * Run:   node lib/commonjs/cli/mvs/mvs-validate examples/mvs/1cbs.mvsj
 */

import { ArgumentParser } from 'argparse';
import fs from 'fs';

import { setFSModule } from '../../mol-util/data-source';
import { MVSData } from '../../extensions/mvs/mvs-data';


setFSModule(fs);

/** Command line argument values for `main` */
interface Args {
    input: string[],
    no_extra: boolean,
}

/** Return parsed command line arguments for `main` */
function parseArguments(): Args {
    const parser = new ArgumentParser({ description: 'Command-line application for validating MolViewSpec files. Prints validation status (OK/FAILED) to stdout, detailed validation issues to stderr. Exits with a zero exit code if all input files are OK.' });
    parser.add_argument('input', { nargs: '+', help: 'Input file(s) in .mvsj format' });
    parser.add_argument('--no-extra', { action: 'store_true', help: 'Treat presence of extra node params as an issue.' });
    const args = parser.parse_args();
    return { ...args };
}

/** Main workflow for validating MolViewSpec files. Returns the number of failed input files. */
function main(args: Args): number {
    let nFailed = 0;
    for (const input of args.input) {
        const data = fs.readFileSync(input, { encoding: 'utf8' });
        const mvsData = MVSData.fromMVSJ(data);
        const issues = MVSData.validationIssues(mvsData, { noExtra: args.no_extra });
        const status = issues ? 'FAILED' : 'OK';
        console.log(`${status.padEnd(6)} ${input}`);
        if (issues) {
            nFailed++;
            for (const issue of issues) {
                console.error(issue);
            }
        }
    }
    return nFailed;
}

const nFailed = main(parseArguments());
if (nFailed > 0) {
    process.exitCode = 1;
}
