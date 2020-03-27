/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';
import * as path from 'path';
import { JobManager, Job, JobEntry } from './jobs';
import { ConsoleLogger } from '../../../mol-util/console-logger';
import { resolveJob } from './query';
import { StructureCache } from './structure-wrapper';
import { now } from '../../../mol-util/now';
import { PerformanceMonitor } from '../../../mol-util/performance-monitor';
import { QueryName } from './api';
import { makeDir } from '../../../mol-util/make-dir';

export type LocalInput = {
    input: string,
    output: string,
    query: QueryName,
    modelNums?: number[],
    params?: any,
    binary?: boolean
}[];

export async function runLocal(input: LocalInput) {
    if (!input.length) {
        ConsoleLogger.error('Local', 'No input');
        return;
    }

    for (const job of input) {
        const binary = /\.bcif/.test(job.output);
        JobManager.add({
            entries: [JobEntry({
                entryId: job.input,
                queryName: job.query,
                queryParams: job.params || { },
                modelNums: job.modelNums,
            })],
            options: {
                outputFilename: job.output,
                binary
            }
        });
    }
    JobManager.sort();

    const started = now();

    let job: Job | undefined = JobManager.getNext();
    let key = job.entries[0].key;
    let progress = 0;
    while (job) {
        try {
            const encoder = await resolveJob(job);
            const writer = wrapFileToWriter(job.outputFilename!);
            encoder.writeTo(writer);
            writer.end();
            ConsoleLogger.logId(job.id, 'Query', 'Written.');

            if (JobManager.hasNext()) {
                job = JobManager.getNext();
                if (key !== job.entries[0].key) StructureCache.expire(key);
                key = job.entries[0].key;
            } else {
                break;
            }
        } catch (e) {
            ConsoleLogger.errorId(job.id, e);
        }
        ConsoleLogger.log('Progress', `[${++progress}/${input.length}] after ${PerformanceMonitor.format(now() - started)}.`);
    }

    ConsoleLogger.log('Progress', `Done in ${PerformanceMonitor.format(now() - started)}.`);
    StructureCache.expireAll();
}

export function wrapFileToWriter(fn: string) {
    const w = {
        open(this: any) {
            if (this.opened) return;
            makeDir(path.dirname(fn));
            this.file = fs.openSync(fn, 'w');
            this.opened = true;
        },
        writeBinary(this: any, data: Uint8Array) {
            this.open();
            fs.writeSync(this.file, Buffer.from(data.buffer));
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