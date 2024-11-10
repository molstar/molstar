/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';
import * as path from 'path';
import * as express from 'express';
import { ModelServerConfig as Config, ModelServerConfig, mapSourceAndIdToFilename } from '../config';
import { ConsoleLogger } from '../../../mol-util/console-logger';
import { resolveJob } from './query';
import { JobManager, JobEntry, ResultWriterParams } from './jobs';
import { UUID } from '../../../mol-util';
import { QueryDefinition, normalizeRestQueryParams, normalizeRestCommonParams, QueryList } from './api';
import { getApiSchema, shortcutIconLink } from './api-schema';
import { swaggerUiAssetsHandler, swaggerUiIndexHandler } from '../../common/swagger-ui';
import { MultipleQuerySpec, getMultiQuerySpecFilename } from './api-web-multiple';
import { SimpleResponseResultWriter, WebResutlWriter, TarballResponseResultWriter } from '../utils/writer';
import { splitCamelCase } from '../../../mol-util/string';
import { healthCheck } from '../../common/util';

function makePath(p: string) {
    return Config.apiPrefix + '/' + p;
}

const responseMap = new Map<UUID, express.Response>();

async function processNextJob() {
    if (!JobManager.hasNext()) return;

    const job = JobManager.getNext();
    responseMap.delete(job.id);
    const writer = job.writer as WebResutlWriter;

    try {
        await resolveJob(job);
    } catch (e) {
        ConsoleLogger.errorId(job.id, '' + e);
        writer.doError(404, '' + e);
    } finally {
        writer.end();
        ConsoleLogger.logId(job.id, 'Query', 'Finished.');
        setImmediate(processNextJob);
    }
}

export function createResultWriter(response: express.Response, params: ResultWriterParams) {
    const filenameBase = params.entryId && params.queryName
        ? `${params.entryId}_${splitCamelCase(params.queryName.replace(/\s/g, '_'), '-').toLowerCase()}`
        : `result`;
    return new SimpleResponseResultWriter(params.filename || `${filenameBase}.${params.encoding}`, response, params.encoding === 'bcif', params.download);
}

function mapQuery(app: express.Express, queryName: string, queryDefinition: QueryDefinition) {
    function createJob(queryParams: any, req: express.Request, res: express.Response) {
        const entryId = req.params.id;
        const commonParams = normalizeRestCommonParams(req.query);
        const resultWriterParams = { encoding: commonParams.encoding!, download: !!commonParams.download, filename: commonParams.filename!, entryId, queryName };
        const jobId = JobManager.add({
            entries: [JobEntry({
                sourceId: commonParams.data_source || ModelServerConfig.defaultSource,
                entryId,
                queryName: queryName as any,
                queryParams,
                modelNums: commonParams.model_nums,
                copyAllCategories: !!commonParams.copy_all_categories,
                transform: commonParams.transform
            })],
            writer: createResultWriter(res, resultWriterParams),
            options: { binary: commonParams.encoding === 'bcif', encoding: commonParams.encoding }
        });
        responseMap.set(jobId, res);
        if (JobManager.size === 1) processNextJob();
    }

    app.get(makePath('v1/:id/' + queryName), (req, res) => {
        const queryParams = normalizeRestQueryParams(queryDefinition, req.query);
        createJob(queryParams, req, res);
    });

    app.post(makePath('v1/:id/' + queryName), (req, res) => {
        const queryParams = req.body;
        createJob(queryParams, req, res);
    });
}

function serveStatic(req: express.Request, res: express.Response) {
    const source = req.params.source === 'bcif'
        ? 'pdb-bcif'
        : req.params.source === 'cif'
            ? 'pdb-cif'
            : req.params.source;

    const id = req.params.id;
    const [fn, format] = mapSourceAndIdToFilename(source, id);
    const binary = format === 'bcif' || fn.indexOf('.bcif') > 0;

    if (!fn || !fs.existsSync(fn)) {
        res.status(404);
        res.end();
        return;
    }
    fs.readFile(fn, (err, data) => {
        if (err) {
            res.status(404);
            res.end();
            return;
        }

        const f = path.parse(fn);
        res.writeHead(200, {
            'Content-Type': binary ? 'application/octet-stream' : 'text/plain; charset=utf-8',
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Headers': 'X-Requested-With',
            'Content-Disposition': `inline; filename="${f.name}${f.ext}"`
        });
        res.write(data);
        res.end();
    });
}

function createMultiJob(spec: MultipleQuerySpec, res: express.Response) {
    const writer = spec.asTarGz
        ? new TarballResponseResultWriter(getMultiQuerySpecFilename(), res)
        : createResultWriter(res, { encoding: spec.encoding!, download: !!spec.download, filename: spec.filename });

    if (spec.queries.length > ModelServerConfig.maxQueryManyQueries) {
        writer.doError(400, `query-many queries limit (${ModelServerConfig.maxQueryManyQueries}) exceeded.`);
        return;
    }

    const jobId = JobManager.add({
        entries: spec.queries.map(q => JobEntry({
            sourceId: q.data_source || ModelServerConfig.defaultSource,
            entryId: q.entryId,
            queryName: q.query,
            queryParams: q.params || { },
            modelNums: q.model_nums,
            copyAllCategories: !!q.copy_all_categories
        })),
        writer,
        options: { binary: spec.encoding?.toLowerCase() === 'bcif', tarball: spec.asTarGz }
    });
    responseMap.set(jobId, res);
    if (JobManager.size === 1) processNextJob();
}

export function initWebApi(app: express.Express) {
    app.use(express.json({ limit: '1mb' }));

    app.get(makePath('static/:source/:id'), (req, res) => serveStatic(req, res));
    app.get(makePath('v1/static/:source/:id'), (req, res) => serveStatic(req, res));

    app.get(makePath('v1/query-many'), (req, res) => {
        const query = /\?query=(.*)$/.exec(req.url)![1];
        const params = JSON.parse(decodeURIComponent(query));
        createMultiJob(params, res);
    });
    app.post(makePath('v1/query-many'), (req, res) => {
        const params = req.body;
        req.setTimeout(ModelServerConfig.requestTimeoutMs);
        createMultiJob(params, res);
    });

    app.use(express.json({ limit: '20mb' }));

    for (const q of QueryList) {
        mapQuery(app, q.name, q.definition);
    }

    // Reports server health depending on `healthCheckPath` config prop
    app.get(makePath('health-check'), (_, res) => healthCheck(res, ModelServerConfig.healthCheckPath));

    const schema = getApiSchema();

    app.get(makePath('openapi.json'), (req, res) => {
        res.writeHead(200, {
            'Content-Type': 'application/json; charset=utf-8',
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Headers': 'X-Requested-With'
        });
        res.end(JSON.stringify(schema));
    });

    app.use(makePath(''), swaggerUiAssetsHandler());
    app.get(makePath(''), swaggerUiIndexHandler({
        openapiJsonUrl: makePath('openapi.json'),
        apiPrefix: Config.apiPrefix,
        title: 'ModelServer API',
        shortcutIconLink
    }));
}