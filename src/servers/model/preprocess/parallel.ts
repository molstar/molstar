/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as path from 'path'
import * as cluster from 'cluster'
import { now } from 'mol-task';
import { PerformanceMonitor } from 'mol-util/performance-monitor';
import { preprocessFile } from './preprocess';

export interface PreprocessEntry {
    source: string,
    cif?: string,
    bcif?: string
}

export interface ParallelPreprocessConfig {
    numProcesses?: number,
    entries: PreprocessEntry[]
}

export function runMaster(config: ParallelPreprocessConfig) {
    const parts = partitionArray(config.entries, config.numProcesses || 1);
    // const numForks = Math.min(parts.length, config.numProcesses);

    const started = now();
    let progress = 0;
    const onMessage = (msg: any) => {
        if (msg.type === 'tick') {
            progress++;
            const elapsed = now() - started;
            console.log(`[${progress}/${config.entries.length}] in ${PerformanceMonitor.format(elapsed)} (avg ${PerformanceMonitor.format(elapsed / progress)}).`);
        } else if (msg.type === 'error') {
            console.error(`${msg.id}: ${msg.error}`)
        }
    }

    for (const _ of parts) {
        const worker = cluster.fork();
        worker.on('message', onMessage);
    }

    let i = 0;
    for (const id in cluster.workers) {
        cluster.workers[id]!.send(parts[i++]);
    }
}

export function runChild() {
    process.on('message', async (entries: PreprocessEntry[]) => {
        for (const entry of entries) {
            try {
                await preprocessFile(entry.source, entry.cif, entry.bcif);
            } catch (e) {
                process.send!({ type: 'error', id: path.parse(entry.source).name, error: '' + e });
            }
            process.send!({ type: 'tick' });
        }
        process.exit();
    });
}

function partitionArray<T>(xs: T[], count: number): T[][] {
    const ret: T[][] = [];
    const s = Math.ceil(xs.length / count);
    for (let i = 0; i < xs.length; i += s) {
        const bucket: T[] = [];
        for (let j = i, _j = Math.min(xs.length, i + s); j < _j; j++) {
            bucket.push(xs[j]);
        }
        ret.push(bucket);
    }
    return ret;
}

