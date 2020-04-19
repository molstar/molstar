/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from '../../../mol-util';
import { getQueryByName, QueryDefinition, QueryName, QueryParams } from './api';
import { LinkedList } from '../../../mol-data/generic';
import { ResultWriter } from '../utils/writer';

export interface ResponseFormat {
    tarball: boolean,
    isBinary: boolean
}

export interface Job {
    id: UUID,
    datetime_utc: string,

    entries: JobEntry[],

    responseFormat: ResponseFormat,
    outputFilename?: string,

    writer: ResultWriter
}

export interface JobDefinition {
    entries: JobEntry[],
    writer: ResultWriter,
    options?: { outputFilename?: string, binary?: boolean, tarball?: boolean }
}

export interface JobEntry {
    job: Job,
    sourceId: '_local_' | string,
    entryId: string,
    key: string,

    queryDefinition: QueryDefinition,
    normalizedParams: any,
    modelNums?: number[],
    copyAllCategories: boolean
}

interface JobEntryDefinition<Name extends QueryName> {
    sourceId?: string, // = '_local_',
    entryId: string,
    queryName: Name,
    queryParams: QueryParams<Name>,
    modelNums?: number[],
    copyAllCategories: boolean
}

export function JobEntry<Name extends QueryName>(definition: JobEntryDefinition<Name>): JobEntry {
    const queryDefinition = getQueryByName(definition.queryName);
    if (!queryDefinition) throw new Error(`Query '${definition.queryName}' is not supported.`);

    const normalizedParams = definition.queryParams;
    const sourceId = definition.sourceId || '_local_';

    return {
        job: void 0 as any,
        key: `${sourceId}/${definition.entryId}`,
        sourceId,
        entryId: definition.entryId,
        queryDefinition,
        normalizedParams,
        modelNums: definition.modelNums,
        copyAllCategories: !!definition.copyAllCategories
    };
}

export function createJob(definition: JobDefinition): Job {
    const job: Job = {
        id: UUID.create22(),
        datetime_utc: `${new Date().toISOString().replace(/T/, ' ').replace(/\..+/, '')}`,
        entries: definition.entries,
        writer: definition.writer,
        responseFormat: { isBinary: !!(definition.options && definition.options.binary), tarball: !!definition?.options?.tarball },
        outputFilename: definition.options && definition.options.outputFilename
    };
    definition.entries.forEach(e => e.job = job);
    return job;
}

class _JobQueue {
    private list: LinkedList<Job> = LinkedList();

    get size() {
        return this.list.count;
    }

    add(definition: JobDefinition) {
        const job = createJob(definition);
        this.list.addLast(job);
        return job.id;
    }

    hasNext(): boolean {
        return this.list.count > 0;
    }

    getNext(): Job {
        return this.list.removeFirst()!;
    }

    /** Sort the job list by key = sourceId/entryId */
    sort() {
        if (this.list.count === 0) return;

        const jobs: Job[] = [];
        for (let j = this.list.first; !!j; j = j.next) {
            jobs[jobs.length] = j.value;
        }

        jobs.sort((a, b) => a.entries[0]?.key < b.entries[0]?.key ? -1 : 1);

        this.list = LinkedList();
        for (const j of jobs) {
            this.list.addLast(j);
        }
    }
}

export const JobManager = new _JobQueue();