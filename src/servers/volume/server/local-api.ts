/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Api from './api';
import * as Data from './query/data-model';
import * as Coordinate from './algebra/coordinate';

import * as fs from 'fs';
import * as path from 'path';

export interface JobEntry {
    source: {
        filename: string,
        name: string,
        id: string
    },
    query: {
        kind: 'box' | 'cell',
        space?: 'fractional' | 'cartesian',
        bottomLeft?: number[],
        topRight?: number[],
    }
    params: {
        /** Determines the detail level as specified in server-config */
        detail?: number,
        /**
         * Determines the sampling level:
         * 1: Original data
         * 2: Downsampled by factor 1/2
         * ...
         * N: downsampled 1/2^(N-1)
         */
        forcedSamplingLevel?: number,
        asBinary: boolean,
    },
    outputFolder: string
}

export async function run(jobs: JobEntry[]) {
    let progress = 0;
    let started = getTime();
    for (const job of jobs) {
        try {
            await query(job);
        } catch (e) {
            console.error(e);
        }
        progress++;
        const elapsed = (getTime() - started) / 1000;
        console.log(`[Progress] ${progress}/${jobs.length} in ${elapsed.toFixed(2)}s`);
    }
}

function getTime() {
    let t = process.hrtime();
    return t[0] * 1000 + t[1] / 1000000;
}

async function query(job: JobEntry) {
    let box: Data.QueryParamsBox;

    if (job.query.kind.toLocaleLowerCase() === 'cell') {
        box = { kind: 'Cell' };
    } else if (job.query.space === 'fractional') {
        box = {
            kind: 'Fractional',
            a: Coordinate.fractional(job.query.bottomLeft![0], job.query.bottomLeft![1], job.query.bottomLeft![2]),
            b: Coordinate.fractional(job.query.topRight![0], job.query.topRight![1], job.query.topRight![2]),
        };
    } else {
        box = {
            kind: 'Cartesian',
            a: Coordinate.cartesian(job.query.bottomLeft![0], job.query.bottomLeft![1], job.query.bottomLeft![2]),
            b: Coordinate.cartesian(job.query.topRight![0], job.query.topRight![1], job.query.topRight![2]),
        };
    }

    const params: Data.QueryParams = {
        sourceFilename: job.source.filename,
        sourceId: job.source.id,
        asBinary: job.params.asBinary,
        box,
        detail: !job.params.detail ? 0 : job.params.detail,
        forcedSamplingLevel: job.params.forcedSamplingLevel
    };

    if (!fs.existsSync(job.outputFolder)) {
        makeDir(job.outputFolder);
    }

    const filename = path.join(job.outputFolder, Api.getOutputFilename(job.source.name, job.source.id, params));
    const res = () => wrapFile(filename);
    await Api.queryBox(params, res);
}

function makeDir(path: string, root?: string): boolean {
    let dirs = path.split(/\/|\\/g),
        dir = dirs.shift();

    root = (root || '') + dir + '/';

    try {
        fs.mkdirSync(root);
    } catch (e) {
        if (!fs.statSync(root).isDirectory()) throw new Error(e);
    }

    return !dirs.length || makeDir(dirs.join('/'), root);
}

function wrapFile(fn: string) {
    const w = {
        open(this: any) {
            if (this.opened) return;
            this.file = fs.openSync(fn, 'w');
            this.opened = true;
        },
        writeBinary(this: any, data: Uint8Array) {
            this.open();
            fs.writeSync(this.file, Buffer.from(data));
            return true;
        },
        writeString(this: any, data: string) {
            this.open();
            fs.writeSync(this.file, data);
            return true;
        },
        end(this: any) {
            if (!this.opened || this.ended) return;
            fs.close(this.file, function () { });
            this.ended = true;
        },
        file: 0,
        ended: false,
        opened: false
    };

    return w;
}