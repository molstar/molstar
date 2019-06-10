/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as argparse from 'argparse'
import * as LocalApi from './server/local-api'
import VERSION from './server/version'
import * as fs from 'fs'
import { LimitsConfig, addLimitsArgs, setLimitsConfig } from './config';

console.log(`VolumeServer Local ${VERSION}, (c) 2018-2019, Mol* contributors`);
console.log();

function description() {
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

    return `Usage: node local jobs.json\n\nExample jobs.json: ${JSON.stringify(exampleJobs, null, 2)}`
}

const parser = new argparse.ArgumentParser({
    addHelp: true,
    description: description()
});
addLimitsArgs(parser)
parser.addArgument(['jobs'], {
    help: `Path to jobs JSON file.`
})

const config: LimitsConfig & { jobs: string } = parser.parseArgs()
setLimitsConfig(config) // sets the config for global use

async function run() {
    let jobs: LocalApi.JobEntry[];
    try {
        jobs = JSON.parse(fs.readFileSync(config.jobs, 'utf-8'));
    } catch (e) {
        console.log('Error:');
        console.error(e);
        return;
    }

    await LocalApi.run(jobs);
}

run();