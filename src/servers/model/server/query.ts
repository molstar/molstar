/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from 'mol-util';
import { getQueryByName, normalizeQueryParams, QueryDefinition } from './api';
import { getStructure } from './structure-wrapper';
import Config from '../config';
import { Progress, now } from 'mol-task';
import { ConsoleLogger } from 'mol-util/console-logger';
import Writer from 'mol-io/writer/writer';
import * as Encoder from 'mol-io/writer/cif'
import { encode_mmCIF_categories } from 'mol-model/structure/export/mmcif';
import { Selection } from 'mol-model/structure';
import Version from '../version'

export interface ResponseFormat {
    isBinary: boolean
}

export interface Request {
    id: UUID,

    sourceId: '_local_' | string,
    entryId: string,

    queryDefinition: QueryDefinition,
    normalizedParams: any,
    responseFormat: ResponseFormat
}

export function createRequest(sourceId: '_local_' | string, entryId: string, queryName: string, params: any): Request {
    const queryDefinition = getQueryByName(queryName);
    if (!queryDefinition) throw new Error(`Query '${queryName}' is not supported.`);

    const normalizedParams = normalizeQueryParams(queryDefinition, params);

    return {
        id: UUID.create(),
        sourceId,
        entryId,
        queryDefinition,
        normalizedParams,
        responseFormat: { isBinary: !!params.binary }
    };
}

export async function resolveRequest(req: Request, writer: Writer) {
    ConsoleLogger.logId(req.id, 'Query', 'Starting.');

    const wrappedStructure = await getStructure(req.sourceId, req.entryId);
    const structure = req.queryDefinition.structureTransform
        ? await req.queryDefinition.structureTransform(req.normalizedParams, wrappedStructure.structure)
        : wrappedStructure.structure;
    const query = req.queryDefinition.query(req.normalizedParams, structure);
    const result = Selection.unionStructure(await query(structure).run(abortingObserver, 250));

    ConsoleLogger.logId(req.id, 'Query', 'Query finished.');

    const encoder = Encoder.create({ binary: req.responseFormat.isBinary, encoderName: `ModelServer ${Version}` });
    encoder.startDataBlock('result');
    encode_mmCIF_categories(encoder, result);

    ConsoleLogger.logId(req.id, 'Query', 'Encoded.');

    encoder.writeTo(writer);

    ConsoleLogger.logId(req.id, 'Query', 'Written.');
}

const maxTime = Config.maxQueryTimeInMs;
export function abortingObserver(p: Progress) {
    if (now() - p.root.progress.startedTime > maxTime) {
        p.requestAbort(`Exceeded maximum allowed time for a query (${maxTime}ms)`);
    }
}