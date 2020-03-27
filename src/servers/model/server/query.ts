/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../../mol-data/db';
import { CifWriter } from '../../../mol-io/writer/cif';
import { StructureQuery, StructureSelection, Structure } from '../../../mol-model/structure';
import { encode_mmCIF_categories } from '../../../mol-model/structure/export/mmcif';
import { Progress } from '../../../mol-task';
import { now } from '../../../mol-util/now';
import { ConsoleLogger } from '../../../mol-util/console-logger';
import { PerformanceMonitor } from '../../../mol-util/performance-monitor';
import { ModelServerConfig as Config } from '../config';
import Version from '../version';
import { Job, JobEntry } from './jobs';
import { createStructureWrapperFromJobEntry, StructureWrapper, resolveStructures } from './structure-wrapper';
import CifField = CifWriter.Field
import { createModelPropertiesProviderFromConfig, ModelPropertiesProvider } from '../property-provider';

export interface Stats {
    structure: StructureWrapper,
    queryTimeMs: number,
    encodeTimeMs: number,
    resultSize: number
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

    const encoder = CifWriter.createEncoder({
        binary: job.responseFormat.isBinary,
        encoderName: `ModelServer ${Version}`,
        binaryAutoClassifyEncoding: true
    });

    // TODO: how to handle missing entries?
    for (const entry of job.entries) {
        const structure = await createStructureWrapperFromJobEntry(entry, propertyProvider());

        // TODO: this should be unique in case the same structure is queried twice
        // const data = (entry.sourceId === '_local_' ? path.basename(entry.entryId) : entry.entryId).replace(/[^a-z0-9\_]/ig, '').toUpperCase();
        encoder.startDataBlock(structure.cifFrame.header);
        await resolveJobEntry(entry, structure, encoder);
    }

    ConsoleLogger.logId(job.id, 'Query', 'Encoding.');
    encoder.encode();

    return encoder;
}

async function resolveJobEntry(entry: JobEntry, structure: StructureWrapper, encoder: CifWriter.Encoder<any>) {
    ConsoleLogger.logId(entry.job.id, 'Query', `Start ${entry.key}/${entry.queryDefinition.name}.`);

    try {
        perf.start('query');
        const sourceStructures = await resolveStructures(structure, entry.modelNums);
        if (!sourceStructures.length) throw new Error('Model not available');

        let structures: Structure[] = sourceStructures;

        if (entry.queryDefinition.structureTransform) {
            structures = [];
            for (const s of sourceStructures) {
                structures.push(await entry.queryDefinition.structureTransform(entry.normalizedParams, s));
            }
        }

        const queries = structures.map(s => entry.queryDefinition.query(entry.normalizedParams, s));
        const result: Structure[] = [];
        for (let i = 0; i < structures.length; i++) {
            const s = await StructureSelection.unionStructure(StructureQuery.run(queries[i], structures[i], { timeoutMs: Config.queryTimeoutMs }))
            if (s.elementCount > 0) result.push(s);
        }
        perf.end('query');

        ConsoleLogger.logId(entry.job.id, 'Query', `Queried ${entry.key}/${entry.queryDefinition.name}.`);

        perf.start('encode');

        encoder.binaryEncodingProvider = getEncodingProvider(structure);

        // TODO: this actually needs to "reversible" in case of error.
        encoder.writeCategory(_model_server_result, entry);
        encoder.writeCategory(_model_server_params, entry);

        if (entry.queryDefinition.filter) encoder.setFilter(entry.queryDefinition.filter);
        if (result.length > 0) encode_mmCIF_categories(encoder, result);
        if (entry.queryDefinition.filter) encoder.setFilter();
        perf.end('encode');

        const stats: Stats = {
            structure: structure,
            queryTimeMs: perf.time('query'),
            encodeTimeMs: perf.time('encode'),
            resultSize: result.reduce((n, s) => n + s.elementCount, 0)
        };

        encoder.writeCategory(_model_server_stats, stats);
        ConsoleLogger.logId(entry.job.id, 'Query', `Written ${entry.key}/${entry.queryDefinition.name}.`);
        return encoder;
    } catch (e) {
        ConsoleLogger.errorId(entry.job.id, e);
        doError(entry, encoder, e);
    } finally {
        encoder.binaryEncodingProvider = void 0;
    }
}

function getEncodingProvider(structure: StructureWrapper) {
    if (!structure.isBinary) return void 0;
    return CifWriter.createEncodingProviderFromCifFrame(structure.cifFrame);
}

function doError(entry: JobEntry, encoder: CifWriter.Encoder<any>, e: any) {
    encoder.writeCategory(_model_server_result, entry);
    encoder.writeCategory(_model_server_params, entry);
    encoder.writeCategory(_model_server_error, '' + e);
}

const maxTime = Config.queryTimeoutMs;
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

const _model_server_result_fields: CifField<any, JobEntry>[] = [
    string<JobEntry>('job_id', ctx => '' + ctx.job.id),
    string<JobEntry>('datetime_utc', ctx => ctx.job.datetime_utc),
    string<JobEntry>('server_version', ctx => Version),
    string<JobEntry>('query_name', ctx => ctx.queryDefinition.name),
    string<JobEntry>('source_id', ctx => ctx.sourceId),
    string<JobEntry>('entry_id', ctx => ctx.entryId),
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
    int32<Stats>('encode_time_ms', ctx => ctx.encodeTimeMs | 0),
    int32<Stats>('element_count', ctx => ctx.resultSize | 0),
];

const _model_server_result: CifWriter.Category<JobEntry> = {
    name: 'model_server_result',
    instance: (job) => CifWriter.categoryInstance(_model_server_result_fields, { data: job, rowCount: 1 })
};

const _model_server_error: CifWriter.Category<string> = {
    name: 'model_server_error',
    instance: (message) => CifWriter.categoryInstance(_model_server_error_fields, { data: message, rowCount: 1 })
};

const _model_server_params: CifWriter.Category<JobEntry> = {
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