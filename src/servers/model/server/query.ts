/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-data/db';
import { CifWriter } from 'mol-io/writer/cif';
import { StructureQuery, StructureSelection, Structure } from 'mol-model/structure';
import { encode_mmCIF_categories } from 'mol-model/structure/export/mmcif';
import { Progress } from 'mol-task';
import { now } from 'mol-util/now';
import { ConsoleLogger } from 'mol-util/console-logger';
import { PerformanceMonitor } from 'mol-util/performance-monitor';
import Config from '../config';
import Version from '../version';
import { Job } from './jobs';
import { createStructureWrapperFromJob, StructureWrapper, resolveStructures } from './structure-wrapper';
import CifField = CifWriter.Field
import { createModelPropertiesProviderFromConfig, ModelPropertiesProvider } from '../property-provider';

export interface Stats {
    structure: StructureWrapper,
    queryTimeMs: number,
    encodeTimeMs: number
}

const perf = new PerformanceMonitor();

let _propertyProvider: ModelPropertiesProvider;
function propertyProvider() {
    if (_propertyProvider) return _propertyProvider;
    _propertyProvider = createModelPropertiesProviderFromConfig() || (() => []);
    return _propertyProvider;
}

export async function resolveJob(job: Job): Promise<CifWriter.Encoder<any>> {
    ConsoleLogger.logId(job.id, 'Query', 'Starting.');

    const wrappedStructure = await createStructureWrapperFromJob(job, propertyProvider());

    try {
        perf.start('query');
        const sourceStructures = await resolveStructures(wrappedStructure, job.modelNums);
        if (!sourceStructures.length) throw new Error('Model not available');

        let structures: Structure[] = sourceStructures;

        if (job.queryDefinition.structureTransform) {
            structures = [];
            for (const s of sourceStructures) {
                structures.push(await job.queryDefinition.structureTransform(job.normalizedParams, s));
            }
        }

        const queries = structures.map(s => job.queryDefinition.query(job.normalizedParams, s));
        const result: Structure[] = [];
        for (let i = 0; i < structures.length; i++) {
            result.push(await StructureSelection.unionStructure(StructureQuery.run(queries[i], structures[i], Config.maxQueryTimeInMs)));
        }
        perf.end('query');

        const encoder = CifWriter.createEncoder({
            binary: job.responseFormat.isBinary,
            encoderName: `ModelServer ${Version}`,
            binaryEncodingPovider: getEncodingProvider(wrappedStructure),
            binaryAutoClassifyEncoding: true
        });

        ConsoleLogger.logId(job.id, 'Query', 'Query finished.');

        perf.start('encode');
        encoder.startDataBlock(sourceStructures[0].models[0].label.toUpperCase());
        encoder.writeCategory(_model_server_result, job);
        encoder.writeCategory(_model_server_params, job);

        // encoder.setFilter(mmCIF_Export_Filters.onlyPositions);
        encode_mmCIF_categories(encoder, result);
        // encoder.setFilter();
        perf.end('encode');

        const stats: Stats = {
            structure: wrappedStructure,
            queryTimeMs: perf.time('query'),
            encodeTimeMs: perf.time('encode')
        };

        encoder.writeCategory(_model_server_stats, stats);
        encoder.encode();
        ConsoleLogger.logId(job.id, 'Query', 'Encoded.');
        return encoder;
    } catch (e) {
        ConsoleLogger.errorId(job.id, e);
        return doError(job, e);
    }
}

function getEncodingProvider(structure: StructureWrapper) {
    if (!structure.isBinary) return void 0;
    return CifWriter.createEncodingProviderFromCifFrame(structure.cifFrame);
}

function doError(job: Job, e: any) {
    const encoder = CifWriter.createEncoder({ binary: job.responseFormat.isBinary, encoderName: `ModelServer ${Version}` });
    encoder.writeCategory(_model_server_result, job);
    encoder.writeCategory(_model_server_params, job);
    encoder.writeCategory(_model_server_error, '' + e);
    encoder.encode();
    return encoder;
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

const _model_server_result_fields: CifField<any, Job>[] = [
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

const _model_server_error_fields: CifField<number, string>[] = [
    string<string>('message', (ctx, i) => ctx)
];

const _model_server_stats_fields: CifField<number, Stats>[] = [
    int32<Stats>('io_time_ms', ctx => ctx.structure.info.readTime | 0),
    int32<Stats>('parse_time_ms', ctx => ctx.structure.info.parseTime | 0),
    // int32<Stats>('attach_props_time_ms', ctx => ctx.structure.info.attachPropsTime | 0),
    int32<Stats>('create_model_time_ms', ctx => ctx.structure.info.createModelTime | 0),
    int32<Stats>('query_time_ms', ctx => ctx.queryTimeMs | 0),
    int32<Stats>('encode_time_ms', ctx => ctx.encodeTimeMs | 0)
];

const _model_server_result: CifWriter.Category<Job> = {
    name: 'model_server_result',
    instance: (job) => CifWriter.categoryInstance(_model_server_result_fields, { data: job, rowCount: 1 })
};

const _model_server_error: CifWriter.Category<string> = {
    name: 'model_server_error',
    instance: (message) => CifWriter.categoryInstance(_model_server_error_fields, { data: message, rowCount: 1 })
};

const _model_server_params: CifWriter.Category<Job> = {
    name: 'model_server_params',
    instance(job) {
        const params: string[][] = [];
        for (const k of Object.keys(job.normalizedParams)) {
            params.push([k, JSON.stringify(job.normalizedParams[k])]);
        }
        return CifWriter.categoryInstance(_model_server_params_fields, { data: params, rowCount: params.length });
    }
};


const _model_server_stats: CifWriter.Category<Stats> = {
    name: 'model_server_stats',
    instance: (stats) => CifWriter.categoryInstance(_model_server_stats_fields, { data: stats, rowCount: 1 })
}