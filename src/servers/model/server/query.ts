/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from 'mol-util';
import { getQueryByName, normalizeQueryParams, QueryDefinition } from './api';
import { getStructure, StructureWrapper } from './structure-wrapper';
import Config from '../config';
import { Progress, now } from 'mol-task';
import { ConsoleLogger } from 'mol-util/console-logger';
import Writer from 'mol-io/writer/writer';
import { CifWriter } from 'mol-io/writer/cif'
import { encode_mmCIF_categories } from 'mol-model/structure/export/mmcif';
import { Selection } from 'mol-model/structure';
import Version from '../version'
import { Column } from 'mol-data/db';
import { PerformanceMonitor } from 'mol-util/performance-monitor';

export interface ResponseFormat {
    isBinary: boolean
}

export interface Request {
    id: UUID,
    datetime_utc: string,

    sourceId: '_local_' | string,
    entryId: string,

    queryDefinition: QueryDefinition,
    normalizedParams: any,
    responseFormat: ResponseFormat
}

export interface Stats {
    structure: StructureWrapper,
    queryTimeMs: number,
    encodeTimeMs: number
}

export function createRequest(sourceId: '_local_' | string, entryId: string, queryName: string, params: any): Request {
    const queryDefinition = getQueryByName(queryName);
    if (!queryDefinition) throw new Error(`Query '${queryName}' is not supported.`);

    const normalizedParams = normalizeQueryParams(queryDefinition, params);

    return {
        id: UUID.create(),
        datetime_utc: `${new Date().toISOString().replace(/T/, ' ').replace(/\..+/, '')}`,
        sourceId,
        entryId,
        queryDefinition,
        normalizedParams,
        responseFormat: { isBinary: !!params.binary }
    };
}

const perf = new PerformanceMonitor();

export async function resolveRequest(req: Request, writer: Writer) {
    ConsoleLogger.logId(req.id, 'Query', 'Starting.');

    const wrappedStructure = await getStructure(req.sourceId, req.entryId);

    perf.start('query');
    const structure = req.queryDefinition.structureTransform
        ? await req.queryDefinition.structureTransform(req.normalizedParams, wrappedStructure.structure)
        : wrappedStructure.structure;
    const query = req.queryDefinition.query(req.normalizedParams, structure);
    const result = Selection.unionStructure(await query(structure).run(abortingObserver, 250));
    perf.end('query');

    ConsoleLogger.logId(req.id, 'Query', 'Query finished.');

    const encoder = CifWriter.createEncoder({ binary: req.responseFormat.isBinary, encoderName: `ModelServer ${Version}` });

    perf.start('encode');
    encoder.startDataBlock(structure.units[0].model.label.toUpperCase());
    encoder.writeCategory(_model_server_result, [req]);
    encoder.writeCategory(_model_server_params, [req]);

    // encoder.setFilter(mmCIF_Export_Filters.onlyPositions);
    encode_mmCIF_categories(encoder, result);
    // encoder.setFilter();
    perf.end('encode');

    ConsoleLogger.logId(req.id, 'Query', 'Encoded.');

    const stats: Stats = {
        structure: wrappedStructure,
        queryTimeMs: perf.time('query'),
        encodeTimeMs: perf.time('encode')
    };

    encoder.writeCategory(_model_server_stats, [stats]);
    encoder.encode();

    encoder.writeTo(writer);

    ConsoleLogger.logId(req.id, 'Query', 'Written.');
}

const maxTime = Config.maxQueryTimeInMs;
export function abortingObserver(p: Progress) {
    if (now() - p.root.progress.startedTime > maxTime) {
        p.requestAbort(`Exceeded maximum allowed time for a query (${maxTime}ms)`);
    }
}

import CifField = CifWriter.Field

function string<T>(name: string, str: (data: T, i: number) => string, isSpecified?: (data: T) => boolean): CifField<number, T> {
    if (isSpecified) {
        return CifField.str(name, (i, d) => str(d, i), { valueKind: (i, d) => isSpecified(d) ? Column.ValueKind.Present : Column.ValueKind.NotPresent });
    }
    return CifField.str(name, (i, d) => str(d, i));
}

function int32<T>(name: string, value: (data: T) => number): CifField<number, T> {
    return CifField.int(name, (i, d) => value(d));
}

const _model_server_result_fields: CifField<number, Request>[] = [
    string<Request>('request_id', ctx => '' + ctx.id),
    string<Request>('datetime_utc', ctx => ctx.datetime_utc),
    string<Request>('server_version', ctx => Version),
    string<Request>('query_name', ctx => ctx.queryDefinition.name),
    string<Request>('source_id', ctx => ctx.sourceId),
    string<Request>('entry_id', ctx => ctx.entryId),
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


const _model_server_result: CifWriter.Category<Request> = {
    name: 'model_server_result',
    instance: (request) => ({ data: request, fields: _model_server_result_fields, rowCount: 1 })
};

const _model_server_params: CifWriter.Category<Request> = {
    name: 'model_server_params',
    instance(request) {
        const params: string[][] = [];
        for (const k of Object.keys(request.normalizedParams)) {
            params.push([k, '' + request.normalizedParams[k]]);
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