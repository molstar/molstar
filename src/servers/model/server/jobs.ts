/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from '../../../mol-util';
import { getQueryByName, normalizeQueryParams, QueryDefinition, QueryName, QueryParams } from './api';
import { LinkedList } from '../../../mol-data/generic';

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
    modelNums?: number[],

    outputFilename?: string
}

export interface JobDefinition<Name extends QueryName> {
    sourceId?: string, // = '_local_',
    entryId: string,
    queryName: Name,
    queryParams: QueryParams<Name>,
    options?: { modelNums?: number[], outputFilename?: string, binary?: boolean }
}

export function createJob<Name extends QueryName>(definition: JobDefinition<Name>): Job {
    const queryDefinition = getQueryByName(definition.queryName);
    if (!queryDefinition) throw new Error(`Query '${definition.queryName}' is not supported.`);

    const normalizedParams = normalizeQueryParams(queryDefinition, definition.queryParams);
    const sourceId = definition.sourceId || '_local_';
    return {
        id: UUID.create22(),
        datetime_utc: `${new Date().toISOString().replace(/T/, ' ').replace(/\..+/, '')}`,
        key: `${sourceId}/${definition.entryId}`,
        sourceId,
        entryId: definition.entryId,
        queryDefinition,
        normalizedParams,
        responseFormat: { isBinary: !!(definition.options && definition.options.binary) },
        modelNums: definition.options && definition.options.modelNums,
        outputFilename: definition.options && definition.options.outputFilename
    };
}

class _JobQueue {
    private list: LinkedList<Job> = LinkedList();

    get size() {
        return this.list.count;
    }

    add<Name extends QueryName>(definition: JobDefinition<Name>) {
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

        jobs.sort((a, b) => a.key < b.key ? -1 : 1);

        this.list = LinkedList();
        for (const j of jobs) {
            this.list.addLast(j);
        }
    }
}

export const JobManager = new _JobQueue();