/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'fs';
import * as path from 'path';
import * as express from 'express';
import Config from '../config';
import { ConsoleLogger } from 'mol-util/console-logger';
import { resolveJob } from './query';
import { JobManager } from './jobs';
import { UUID } from 'mol-util';
import { LandingPage } from './landing';

function makePath(p: string) {
    return Config.appPrefix + '/' + p;
}

function wrapResponse(fn: string, res: express.Response) {
    const w = {
        doError(this: any, code = 404, message = 'Not Found.') {
            if (!this.headerWritten) {
                res.status(code).send(message);
                this.headerWritten = true;
            }
            this.end();
        },
        writeHeader(this: any, binary: boolean) {
            if (this.headerWritten) return;
            res.writeHead(200, {
                'Content-Type': binary ? 'application/octet-stream' : 'text/plain; charset=utf-8',
                'Access-Control-Allow-Origin': '*',
                'Access-Control-Allow-Headers': 'X-Requested-With',
                'Content-Disposition': `inline; filename="${fn}"`
            });
            this.headerWritten = true;
        },
        writeBinary(this: any, data: Uint8Array) {
            if (!this.headerWritten) this.writeHeader(true);
            return res.write(new Buffer(data.buffer));
        },
        writeString(this: any, data: string) {
            if (!this.headerWritten) this.writeHeader(false);
            return res.write(data);
        },
        end(this: any) {
            if (this.ended) return;
            res.end();
            this.ended = true;
        },
        ended: false,
        headerWritten: false
    };

    return w;
}

const responseMap = new Map<UUID, express.Response>();

async function processNextJob() {
    if (!JobManager.hasNext()) return;

    const job = JobManager.getNext();
    const response = responseMap.get(job.id)!;
    responseMap.delete(job.id);

    const filenameBase = `${job.entryId}_${job.queryDefinition.name.replace(/\s/g, '_')}`
    const writer = wrapResponse(job.responseFormat.isBinary ? `${filenameBase}.bcif` : `${filenameBase}.cif`, response);

    try {
        const encoder = await resolveJob(job);
        writer.writeHeader(job.responseFormat.isBinary);
        encoder.writeTo(writer);
    } catch (e) {
        ConsoleLogger.errorId(job.id, '' + e);
        writer.doError(404, '' + e);
    } finally {
        writer.end();
        ConsoleLogger.logId(job.id, 'Query', 'Finished.');
        setImmediate(processNextJob);
    }
}

// function mapQuery(app: express.Express, queryName: string, queryDefinition: QueryDefinition) {
//     app.get(makePath(':entryId/' + queryName), (req, res) => {
//         ConsoleLogger.log('Server', `Query '${req.params.entryId}/${queryName}'...`);

//         if (JobManager.size >= Config.maxQueueLength) {
//             res.status(503).send('Too many queries, please try again later.');
//             res.end();
//             return;
//         }

//         const jobId = JobManager.add('pdb', req.params.entryId, queryName, req.query);
//         responseMap.set(jobId, res);
//         if (JobManager.size === 1) processNextJob();
//     });
// }

export function initWebApi(app: express.Express) {
    app.get(makePath('static/:format/:id'), async (req, res) => {
        const binary = req.params.format === 'bcif';
        const id = req.params.id;
        const fn = Config.mapFile(binary ? 'pdb-bcif' : 'pdb-cif', id);
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
    })

    app.get(makePath('api/v1'), (req, res) => {
        const query = /\?(.*)$/.exec(req.url)![1];
        const args = JSON.parse(decodeURIComponent(query));
        const name = args.name;
        const entryId = args.id;
        const params = args.params || { };
        const jobId = JobManager.add('pdb', entryId, name, params, args.modelNums);
        responseMap.set(jobId, res);
        if (JobManager.size === 1) processNextJob();
    });

    app.get('*', (req, res) => {
        res.send(LandingPage);
    });

    // for (const q of QueryList) {
    //     mapQuery(app, q.name, q.definition);
    // }
}