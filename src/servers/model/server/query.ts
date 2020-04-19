/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as path from 'path';
import { Column } from '../../../mol-data/db';
import { CifWriter } from '../../../mol-io/writer/cif';
import { Structure, StructureQuery, StructureSelection } from '../../../mol-model/structure';
import { encode_mmCIF_categories } from '../../../mol-model/structure/export/mmcif';
import { Progress } from '../../../mol-task';
import { ConsoleLogger } from '../../../mol-util/console-logger';
import { now } from '../../../mol-util/now';
import { PerformanceMonitor } from '../../../mol-util/performance-monitor';
import { ModelServerConfig as Config } from '../config';
import { createModelPropertiesProviderFromConfig, ModelPropertiesProvider } from '../property-provider';
import Version from '../version';
import { Job, JobEntry } from './jobs';
import { createStructureWrapperFromJobEntry, resolveStructures, StructureWrapper } from './structure-wrapper';
import CifField = CifWriter.Field
import { splitCamelCase } from '../../../mol-util/string';

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

export async function resolveJob(job: Job) {
    if (job.responseFormat.tarball) {
        return resolveMultiFile(job);
    } else {
        return resolveSingleFile(job);
    }
}

async function resolveSingleFile(job: Job) {
    ConsoleLogger.logId(job.id, 'Query', 'Starting.');

    const encoder = CifWriter.createEncoder({
        binary: job.responseFormat.isBinary,
        encoderName: `ModelServer ${Version}`,
        binaryAutoClassifyEncoding: true
    });

    const headerMap = new Map<string, number>();

    for (const entry of job.entries) {
        let hasDataBlock = false;
        try {
            const structure = await createStructureWrapperFromJobEntry(entry, propertyProvider());

            let header = structure.cifFrame.header.toUpperCase();
            if (headerMap.has(header)) {
                const i = headerMap.get(header)! + 1;
                headerMap.set(header, i);
                header += ' ' + i;
            } else {
                headerMap.set(header, 0);
            }

            encoder.startDataBlock(header);
            hasDataBlock = true;
            await resolveJobEntry(entry, structure, encoder);
        } catch (e) {
            if (job.entries.length === 1) {
                throw e;
            } else {
                if (!hasDataBlock) {
                    createErrorDataBlock(entry, encoder);
                }
                doError(entry, encoder, e);
                ConsoleLogger.errorId(entry.job.id, '' + e);
            }
        }
    }

    ConsoleLogger.logId(job.id, 'Query', 'Encoding.');
    encoder.encode();
    encoder.writeTo(job.writer);
}

function getFilename(i: number, entry: JobEntry, header: string, isBinary: boolean) {
    return `${i}_${header.toLowerCase()}_${splitCamelCase(entry.queryDefinition.name.replace(/\s/g, '_'), '-').toLowerCase()}.${isBinary ? 'bcif' : 'cif'}`;
}

async function resolveMultiFile(job: Job) {
    ConsoleLogger.logId(job.id, 'Query', 'Starting.');

    let i = 0;
    for (const entry of job.entries) {

        const encoder = CifWriter.createEncoder({
            binary: job.responseFormat.isBinary,
            encoderName: `ModelServer ${Version}`,
            binaryAutoClassifyEncoding: true
        });

        let hasDataBlock = false;
        let header = '';
        try {
            const structure = await createStructureWrapperFromJobEntry(entry, propertyProvider());
            header = structure.cifFrame.header;
            encoder.startDataBlock(structure.cifFrame.header);
            hasDataBlock = true;
            await resolveJobEntry(entry, structure, encoder);
        } catch(e) {
            if (!hasDataBlock) {
                header = createErrorDataBlock(entry, encoder);
            }
            ConsoleLogger.errorId(entry.job.id, '' + e);
            doError(entry, encoder, e);
        }

        ConsoleLogger.logId(job.id, 'Query', `Encoding ${entry.key}/${entry.queryDefinition.name}`);
        encoder.encode();

        job.writer.beginEntry(getFilename(++i, entry, header, job.responseFormat.isBinary), encoder.getSize());
        encoder.writeTo(job.writer);
        job.writer.endEntry();
        ConsoleLogger.logId(job.id, 'Query', `Written ${entry.key}/${entry.queryDefinition.name}`);

        // await fileEntry;
    }
}

function createErrorDataBlock(job: JobEntry, encoder: CifWriter.Encoder<any>) {
    let header;
    if (job.sourceId === '_local_') header = path.basename(job.entryId).replace(/[^a-z0-9\-]/gi, '').toUpperCase();
    else header = job.entryId.replace(/[^a-z0-9\-]/gi, '').toUpperCase();
    encoder.startDataBlock(header);
    return header;
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
            const s = StructureSelection.unionStructure(StructureQuery.run(queries[i], structures[i], { timeoutMs: Config.queryTimeoutMs }));
            if (s.elementCount > 0) result.push(s);
        }
        perf.end('query');

        ConsoleLogger.logId(entry.job.id, 'Query', `Queried ${entry.key}/${entry.queryDefinition.name}.`);

        perf.start('encode');

        encoder.binaryEncodingProvider = getEncodingProvider(structure);

        // TODO: this actually needs to "reversible" in case of error.
        encoder.writeCategory(_model_server_result, entry);
        encoder.writeCategory(_model_server_params, entry);

        if (!entry.copyAllCategories && entry.queryDefinition.filter) encoder.setFilter(entry.queryDefinition.filter);
        if (result.length > 0) encode_mmCIF_categories(encoder, result, { copyAllCategories: entry.copyAllCategories });
        if (!entry.copyAllCategories && entry.queryDefinition.filter) encoder.setFilter();
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
};