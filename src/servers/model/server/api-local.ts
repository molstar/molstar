/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ConsoleLogger } from '../../../mol-util/console-logger';
import { now } from '../../../mol-util/now';
import { PerformanceMonitor } from '../../../mol-util/performance-monitor';
import { FileResultWriter } from '../utils/writer';
import { QueryName } from './api';
import { Job, JobEntry, JobManager } from './jobs';
import { resolveJob } from './query';
import { StructureCache } from './structure-wrapper';

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
            writer: new FileResultWriter(job.output),
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
            const writer = job.writer;
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