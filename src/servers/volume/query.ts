#!/usr/bin/env node
/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'fs';
import { configureLocal } from './config';
import * as LocalApi from './server/local-api';

const config = configureLocal();

if (config.jobsTemplate !== null) {
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
    console.log(JSON.stringify(exampleJobs, null, 2));
    process.exit();
}

async function run() {
    let jobs: LocalApi.JobEntry[];
    try {
        if (!config.jobs) {
            throw new Error(`Please provide 'jobs' argument. See [-h] for help.`);
        }

        jobs = JSON.parse(fs.readFileSync(config.jobs, 'utf-8'));
    } catch (e) {
        console.error('' + e);
        return;
    }

    await LocalApi.run(jobs);
}

run();