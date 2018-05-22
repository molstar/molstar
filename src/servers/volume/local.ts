/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as LocalApi from './server/local-api'
import VERSION from './server/version'

import * as fs from 'fs'

console.log(`VolumeServer ${VERSION}, (c) 2016 - now, David Sehnal`);
console.log();

function help() {
    const exampleJobs: LocalApi.JobEntry[] = [{
        source: {
            filename: `g:/test/mdb/xray-1tqn.mdb`,
            name: 'xray',
            id: '1tqn',
        },
        query: {
            kind: 'box',
            space: 'cartesian',
            bottomLeft: [-42.996, -64.169, -45.335],
            topRight: [8.768, 15.316, 21.599]
        },
        params: {
            forcedSamplingLevel: 2,
            asBinary: true
        },
        outputFolder: 'g:/test/local-test'
    }, {
        source: {
            filename: `g:/test/mdb/emd-8116.mdb`,
            name: 'em',
            id: '8116',
        },
        query: {
            kind: 'cell'
        },
        params: {
            detail: 4,
            asBinary: true
        },
        outputFolder: 'g:/test/local-test'
    }];

    console.log('Usage: node local jobs.json');
    console.log();
    console.log('Example jobs.json:');
    console.log(JSON.stringify(exampleJobs, null, 2));
}

async function run() {
    if (process.argv.length !== 3) {
        help();
        return;
    }

    let jobs: LocalApi.JobEntry[];
    try {
        jobs = JSON.parse(fs.readFileSync(process.argv[2], 'utf-8'));
    } catch (e) {
        console.log('Error:');
        console.error(e);
        return;
    }

    await LocalApi.run(jobs);
}

run();