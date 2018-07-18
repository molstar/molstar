/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-data/db';
import { CifWriter } from 'mol-io/writer/cif';
import Writer from 'mol-io/writer/writer';
import { StructureQuery, StructureSelection } from 'mol-model/structure';
import { encode_mmCIF_categories } from 'mol-model/structure/export/mmcif';
import { now, Progress } from 'mol-task';
import { ConsoleLogger } from 'mol-util/console-logger';
import { PerformanceMonitor } from 'mol-util/performance-monitor';
import Config from '../config';
import Version from '../version';
import { Job } from './jobs';
import { getStructure, StructureWrapper } from './structure-wrapper';
import CifField = CifWriter.Field

export interface Stats {
    structure: StructureWrapper,
    queryTimeMs: number,
    encodeTimeMs: number
}

const perf = new PerformanceMonitor();

export async function resolveJob(job: Job, writer: Writer) {
    ConsoleLogger.logId(job.id, 'Query', 'Starting.');

    const wrappedStructure = await getStructure(job);

    perf.start('query');
    const structure = job.queryDefinition.structureTransform
        ? await job.queryDefinition.structureTransform(job.normalizedParams, wrappedStructure.structure)
        : wrappedStructure.structure;
    const query = job.queryDefinition.query(job.normalizedParams, structure);
    const result = StructureSelection.unionStructure(StructureQuery.run1(query, structure));
    perf.end('query');

    ConsoleLogger.logId(job.id, 'Query', 'Query finished.');

    const encoder = CifWriter.createEncoder({ binary: job.responseFormat.isBinary, encoderName: `ModelServer ${Version}` });

    perf.start('encode');
    encoder.startDataBlock(structure.units[0].model.label.toUpperCase());
    encoder.writeCategory(_model_server_result, [job]);
    encoder.writeCategory(_model_server_params, [job]);

    // encoder.setFilter(mmCIF_Export_Filters.onlyPositions);
    encode_mmCIF_categories(encoder, result);
    // encoder.setFilter();
    perf.end('encode');

    ConsoleLogger.logId(job.id, 'Query', 'Encoded.');

    const stats: Stats = {
        structure: wrappedStructure,
        queryTimeMs: perf.time('query'),
        encodeTimeMs: perf.time('encode')
    };

    encoder.writeCategory(_model_server_stats, [stats]);
    encoder.encode();

    encoder.writeTo(writer);

    ConsoleLogger.logId(job.id, 'Query', 'Written.');
}

const maxTime = Config.maxQueryTimeInMs;
export function abortingObserver(p: Progress) {
    if (now() - p.root.progress.startedTime > maxTime) {
        p.requestAbort(`Exceeded maximum allowed time for a query (${maxTime}ms)`);
    }
}

function string<T>(name: string, str: (data: T, i: number) => string, isSpecified?: (data: T) => boolean): CifField<number, T> {
    if (isSpecified) {
        return CifField.str(name, (i, d) => str(d, i), { valueKind: (i, d) => isSpecified(d) ? Column.ValueKind.Present : Column.ValueKind.NotPresent });
    }
    return CifField.str(name, (i, d) => str(d, i));
}

function int32<T>(name: string, value: (data: T) => number): CifField<number, T> {
    return CifField.int(name, (i, d) => value(d));
}

const _model_server_result_fields: CifField<number, Job>[] = [
    string<Job>('job_id', ctx => '' + ctx.id),
    string<Job>('datetime_utc', ctx => ctx.datetime_utc),
    string<Job>('server_version', ctx => Version),
    string<Job>('query_name', ctx => ctx.queryDefinition.name),
    string<Job>('source_id', ctx => ctx.sourceId),
    string<Job>('entry_id', ctx => ctx.entryId),
];

const _model_server_params_fields: CifField<number, string[]>[] = [
    string<string[]>('name', (ctx, i) => ctx[i][0]),
    string<string[]>('value', (ctx, i) => ctx[i][1])
];

const _model_server_stats_fields: CifField<number, Stats>[] = [
    int32<Stats>('io_time_ms', ctx => ctx.structure.info.readTime | 0),
    int32<Stats>('parse_time_ms', ctx => ctx.structure.info.parseTime | 0),
    int32<Stats>('create_model_time_ms', ctx => ctx.structure.info.createModelTime | 0),
    int32<Stats>('query_time_ms', ctx => ctx.queryTimeMs | 0),
    int32<Stats>('encode_time_ms', ctx => ctx.encodeTimeMs | 0)
];


const _model_server_result: CifWriter.Category<Job> = {
    name: 'model_server_result',
    instance: (job) => ({ data: job, fields: _model_server_result_fields, rowCount: 1 })
};

const _model_server_params: CifWriter.Category<Job> = {
    name: 'model_server_params',
    instance(job) {
        const params: string[][] = [];
        for (const k of Object.keys(job.normalizedParams)) {
            params.push([k, '' + job.normalizedParams[k]]);
        }
        return {
            data: params,

            fields: _model_server_params_fields,
            rowCount: params.length
        }
    }
};


const _model_server_stats: CifWriter.Category<Stats> = {
    name: 'model_server_stats',
    instance: (stats) => ({ data: stats, fields: _model_server_stats_fields, rowCount: 1 })
}