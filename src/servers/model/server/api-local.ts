/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ConsoleLogger } from '../../../mol-util/console-logger';
import { now } from '../../../mol-util/now';
import { PerformanceMonitor } from '../../../mol-util/performance-monitor';
import { FileResultWriter, TarballFileResultWriter } from '../utils/writer';
import { QueryName, QueryParams } from './api';
import { Job, JobEntry, JobManager } from './jobs';
import { resolveJob } from './query';
import { StructureCache } from './structure-wrapper';

export type Entry<Q extends QueryName = QueryName> = {
    input: string,
    query: Q,
    modelNums?: number[],
    copyAllCategories?: boolean,
    params?: QueryParams<Q>,
}

export type LocalInput = {
    queries: Entry[],
    output: string,
    binary?: boolean,
    asTarGz?: boolean,
    gzipLevel?: number
}[];

export async function runLocal(input: LocalInput) {
    if (!input.length) {
        ConsoleLogger.error('Local', 'No input');
        return;
    }

    for (const job of input) {
        const binary = /\.bcif/.test(job.output);
        JobManager.add({
            entries: job.queries.map(q => JobEntry({
                entryId: q.input,
                queryName: q.query,
                queryParams: q.params || { },
                modelNums: q.modelNums,
                copyAllCategories: !!q.copyAllCategories
            })),
            writer: job.asTarGz
                ? new TarballFileResultWriter(job.output, job.gzipLevel)
                : new FileResultWriter(job.output),
            options: {
                outputFilename: job.output,
                binary,
                tarball: job.asTarGz
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
            await resolveJob(job);
            job.writer.end();
            ConsoleLogger.logId(job.id, 'Query', 'Written.');

            if (job.entries.length > 0) StructureCache.expireAll();

            if (JobManager.hasNext()) {
                job = JobManager.getNext();
                if (key !== job.entries[0].key) StructureCache.expire(key);
                key = job.entries[0].key;
            } else {
                break;
            }
        } catch (e) {
            ConsoleLogger.errorId(job.id, e);

            if (JobManager.hasNext()) {
                job = JobManager.getNext();
                if (key !== job.entries[0].key) StructureCache.expire(key);
                key = job.entries[0].key;
            } else {
                break;
            }
        }
        ConsoleLogger.log('Progress', `[${++progress}/${input.length}] after ${PerformanceMonitor.format(now() - started)}.`);
    }

    ConsoleLogger.log('Progress', `Done in ${PerformanceMonitor.format(now() - started)}.`);
    StructureCache.expireAll();
}