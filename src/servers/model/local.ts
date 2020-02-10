/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs'
import Version from './version';
import { LocalInput, runLocal } from './server/api-local';

console.log(`Mol* ModelServer (${Version}), (c) 2018 Mol* authors`);
console.log(``);

let exampleWorkload: LocalInput = [{
    input: 'c:/test/quick/1tqn.cif',
    output: 'c:/test/quick/localapi/1tqn_full.cif',
    query: 'full', // same as defined in Api/Queries
}, {
    input: 'c:/test/quick/1tqn.cif',
    output: 'c:/test/quick/localapi/1tqn_full.bcif',
    query: 'full',
    params: {}
}, {
    input: 'c:/test/quick/1cbs_updated.cif',
    output: 'c:/test/quick/localapi/1cbs_ligint.cif',
    query: 'residueInteraction', // action is case sensitive
    params: { atom_site: { label_comp_id: 'REA' }, radius: 5 }
}, {
    input: 'c:/test/quick/1cbs_updated.cif', // multiple files that are repeated will only be parsed once
    output: 'c:/test/quick/localapi/1cbs_ligint.bcif',
    query: 'residueInteraction',
    params: { atom_site: { label_comp_id: 'REA' } } // parameters are just a JSON version of the query string
}];


if (process.argv.length !== 3) {
    let help = [
        `Usage: `,
        ``,
        `   node local jobs.json`,
        ``,
        `jobs.json is a JSON version of the WebAPI. Query names are case sensitive.`,
        `The jobs are automatically sorted by inputFilenama and the given file is only loaded once.`,
        `All processing errors are sent to stderr.`,
        ``,
        `Jobs example:`,
        ``,
        JSON.stringify(exampleWorkload, null, 2)
    ];

    console.log(help.join('\n'));
} else {
    try {
        const input = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
        runLocal(input);
    } catch (e) {
        console.error(e);
    }
}

// TODO: write utility that splits jobs into multiple chunks?
