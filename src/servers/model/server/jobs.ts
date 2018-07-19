/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from 'mol-util';
import { getQueryByName, normalizeQueryParams, QueryDefinition } from './api';
import { LinkedList } from 'mol-data/generic';

export interface ResponseFormat {
    isBinary: boolean
}

export interface Job {
    id: UUID,
    datetime_utc: string,

    sourceId: '_local_' | string,
    entryId: string,
    key: string,

    queryDefinition: QueryDefinition,
    normalizedParams: any,
    responseFormat: ResponseFormat,

    outputFilename?: string
}

export function createJob(sourceId: '_local_' | string, entryId: string, queryName: string, params: any, outputFilename?: string): Job {
    const queryDefinition = getQueryByName(queryName);
    if (!queryDefinition) throw new Error(`Query '${queryName}' is not supported.`);

    const normalizedParams = normalizeQueryParams(queryDefinition, params);

    return {
        id: UUID.create(),
        datetime_utc: `${new Date().toISOString().replace(/T/, ' ').replace(/\..+/, '')}`,
        key: `${sourceId}/${entryId}`,
        sourceId,
        entryId,
        queryDefinition,
        normalizedParams,
        responseFormat: { isBinary: !!params.binary },
        outputFilename
    };
}

class _JobQueue {
    private list: LinkedList<Job> = LinkedList();

    get size() {
        return this.list.count;
    }

    add(sourceId: '_local_' | string, entryId: string, queryName: string, params: any, outputFilename?: string) {
        const job = createJob(sourceId, entryId, queryName, params, outputFilename);
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

        jobs.sort((a, b) => a.key < b.key ? -1 : 1);

        this.list = LinkedList();
        for (const j of jobs) {
            this.list.addLast(j);
        }
    }
}

export const JobManager = new _JobQueue();