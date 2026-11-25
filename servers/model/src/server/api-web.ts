/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import * as fs from 'fs';
import * as path from 'path';
import * as express from 'express';
import { ModelServerConfig as Config, ModelServerConfig, mapSourceAndIdToFilename } from '@molstar/model-server/config';
import { ConsoleLogger } from '@molstar/core/util/console-logger';
import { resolveJob } from './query.js';
import { JobEntry, JobRunner, createJob, type Job, type JobDefinition, type ResultWriterParams } from './jobs.js';
import { type QueryDefinition, normalizeRestQueryParams, normalizeRestCommonParams, QueryList } from './api.js';
import { getApiSchema, shortcutIconLink } from './api-schema.js';
import { swaggerUiAssetsHandler, swaggerUiIndexHandler } from '@molstar/common-server/swagger-ui';
import { type MultipleQuerySpec, getMultiQuerySpecFilename } from './api-web-multiple.js';
import {
  SimpleResponseResultWriter,
  type WebResutlWriter,
  TarballResponseResultWriter,
} from '@molstar/model-server/utils/writer';
import { splitCamelCase } from '@molstar/core/util/string';
import { healthCheck } from '@molstar/common-server/util';
import { StructureCache } from './structure-wrapper.js';
import { errorStatus } from './errors.js';

function makePath(p: string) {
  return Config.apiPrefix + '/' + p;
}

const MB = 1024 * 1024;

function stats() {
  const { queued, running } = runner.stats;
  const cache = StructureCache.stats;
  const heap = process.memoryUsage().heapUsed;
  return (
    `queued ${queued}, running ${running}, cache ${cache.count} entries ~${Math.round(cache.approximateSize / MB)} MB, ` +
    `heap ${Math.round(heap / MB)} MB`
  );
}

async function processJob(job: Job) {
  const writer = job.writer as WebResutlWriter;
  try {
    await resolveJob(job);
  } catch (e) {
    ConsoleLogger.errorId(job.id, '' + e);
    writer.doError(errorStatus(e), '' + e);
  } finally {
    writer.end();
    ConsoleLogger.logId(job.id, 'Query', `Finished (${stats()}).`);
  }
}

function skipJob(job: Job) {
  job.writer.end();
  ConsoleLogger.logId(job.id, 'Query', 'Skipped: the client is gone.');
}

const runner = new JobRunner(
  processJob,
  () => ({ maxConcurrentJobs: Config.maxConcurrentJobs, maxQueueLength: Config.maxQueueLength }),
  skipJob,
);

/** Stops accepting queries (503) and resolves when the queued and running ones are done. */
export function drainWebApi() {
  return runner.drain();
}

function enqueue(definition: JobDefinition, res: express.Response) {
  const job = createJob(definition);
  if (runner.add(job)) return;

  const writer = job.writer as WebResutlWriter;
  if (runner.isAccepting) {
    ConsoleLogger.logId(job.id, 'Query', `Refused: the queue is full (${stats()}).`);
    res.setHeader('Retry-After', '10');
    writer.doError(503, 'The server is busy, try again later.');
  } else {
    writer.doError(503, 'The server is shutting down.');
  }
}

function badRequest(res: express.Response, e: unknown) {
  res.status(errorStatus(e) === 500 ? 400 : errorStatus(e)).send('' + e);
}

export function createResultWriter(response: express.Response, params: ResultWriterParams) {
  const filenameBase =
    params.entryId && params.queryName
      ? `${params.entryId}_${splitCamelCase(params.queryName.replace(/\s/g, '_'), '-').toLowerCase()}`
      : `result`;
  return new SimpleResponseResultWriter(
    params.filename || `${filenameBase}.${params.encoding}`,
    response,
    params.encoding === 'bcif',
    params.download,
  );
}

// bodies of query-many specs, and of query POSTs (their parameters)
const jsonSpec = express.json({ limit: '1mb' });
const jsonParams = express.json({ limit: '20mb' });

function mapQuery(app: express.Express, queryName: string, queryDefinition: QueryDefinition) {
  function handle(getQueryParams: () => any, req: express.Request, res: express.Response) {
    let definition: JobDefinition;
    try {
      const entryId = Array.isArray(req.params.id) ? req.params.id[0] : req.params.id;
      const commonParams = normalizeRestCommonParams(req.query);
      const resultWriterParams = {
        encoding: commonParams.encoding!,
        download: !!commonParams.download,
        filename: commonParams.filename!,
        entryId,
        queryName,
      };
      definition = {
        entries: [
          JobEntry({
            sourceId: commonParams.data_source || ModelServerConfig.defaultSource,
            entryId,
            queryName: queryName as any,
            queryParams: getQueryParams(),
            modelNums: commonParams.model_nums,
            copyAllCategories: !!commonParams.copy_all_categories,
            transform: commonParams.transform,
          }),
        ],
        writer: createResultWriter(res, resultWriterParams),
        options: { binary: commonParams.encoding === 'bcif', encoding: commonParams.encoding },
      };
    } catch (e) {
      badRequest(res, e);
      return;
    }
    enqueue(definition, res);
  }

  app.get(makePath('v1/:id/' + queryName), (req, res) => {
    handle(() => normalizeRestQueryParams(queryDefinition, req.query), req, res);
  });

  app.post(makePath('v1/:id/' + queryName), jsonParams, (req, res) => {
    handle(() => req.body, req, res);
  });
}

function serveStatic(req: express.Request, res: express.Response) {
  const _source = Array.isArray(req.params.source) ? req.params.source[0] : req.params.source;
  const source = _source === 'bcif' ? 'pdb-bcif' : _source === 'cif' ? 'pdb-cif' : _source;

  const id = Array.isArray(req.params.id) ? req.params.id[0] : req.params.id;
  const [fn, format] = mapSourceAndIdToFilename(source, id) ?? [];

  if (!fn || !fs.existsSync(fn)) {
    res.status(404);
    res.end();
    return;
  }
  const binary = format === 'bcif' || fn.indexOf('.bcif') > 0;
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
      'Content-Disposition': `inline; filename="${f.name}${f.ext}"`,
    });
    res.write(data);
    res.end();
  });
}

function createMultiJob(spec: MultipleQuerySpec, res: express.Response) {
  if (!spec || typeof spec !== 'object' || !Array.isArray(spec.queries)) {
    res.status(400).send(`query-many expects { queries: [...] }.`);
    return;
  }

  const writer = spec.asTarGz
    ? new TarballResponseResultWriter(getMultiQuerySpecFilename(), res)
    : createResultWriter(res, { encoding: spec.encoding!, download: !!spec.download, filename: spec.filename });

  if (spec.queries.length > ModelServerConfig.maxQueryManyQueries) {
    writer.doError(400, `query-many queries limit (${ModelServerConfig.maxQueryManyQueries}) exceeded.`);
    return;
  }

  let entries: JobEntry[];
  try {
    entries = spec.queries.map((q) =>
      JobEntry({
        sourceId: q.data_source || ModelServerConfig.defaultSource,
        entryId: q.entryId,
        queryName: q.query,
        queryParams: q.params || {},
        modelNums: q.model_nums,
        copyAllCategories: !!q.copy_all_categories,
      }),
    );
  } catch (e) {
    writer.doError(400, '' + e);
    return;
  }

  enqueue(
    { entries, writer, options: { binary: spec.encoding?.toLowerCase() === 'bcif', tarball: spec.asTarGz } },
    res,
  );
}

export function initWebApi(app: express.Express) {
  app.get(makePath('static/:source/:id'), (req, res) => serveStatic(req, res));
  app.get(makePath('v1/static/:source/:id'), (req, res) => serveStatic(req, res));

  app.get(makePath('v1/query-many'), (req, res) => {
    // spec is everything after `?query=`, with JSON that may contain &
    const query = /\?query=(.*)$/.exec(req.url);
    if (!query) {
      res.status(400).send('Missing the query parameter.');
      return;
    }
    let params: MultipleQuerySpec;
    try {
      params = JSON.parse(decodeURIComponent(query[1]));
    } catch (e) {
      res.status(400).send('' + e);
      return;
    }
    createMultiJob(params, res);
  });
  app.post(makePath('v1/query-many'), jsonSpec, (req, res) => {
    const params = req.body;
    req.setTimeout(ModelServerConfig.requestTimeoutMs);
    createMultiJob(params, res);
  });

  for (const q of QueryList) {
    mapQuery(app, q.name, q.definition);
  }

  // Reports server health depending on `healthCheckPath` config prop, unhealthy during shutdown
  app.get(makePath('health-check'), (_, res) => {
    if (!runner.isAccepting) {
      res.status(503).send('Shutting down.');
      return;
    }
    healthCheck(res, ModelServerConfig.healthCheckPath);
  });

  const schema = getApiSchema();

  app.get(makePath('openapi.json'), (req, res) => {
    res.writeHead(200, {
      'Content-Type': 'application/json; charset=utf-8',
      'Access-Control-Allow-Origin': '*',
      'Access-Control-Allow-Headers': 'X-Requested-With',
    });
    res.end(JSON.stringify(schema));
  });

  app.use(makePath(''), swaggerUiAssetsHandler());
  app.get(
    makePath(''),
    swaggerUiIndexHandler({
      openapiJsonUrl: makePath('openapi.json'),
      apiPrefix: Config.apiPrefix,
      title: 'ModelServer API',
      shortcutIconLink,
    }),
  );
}
