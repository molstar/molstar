/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as express from 'express';
import Config from '../config';
import { QueryDefinition, QueryList } from './api';
import { ConsoleLogger } from 'mol-util/console-logger';
import { createRequest, resolveRequest } from './query';

function makePath(p: string) {
    return Config.appPrefix + '/' + p;
}

function wrapResponse(fn: string, res: express.Response) {
    const w = {
        do404(this: any) {
            if (!this.headerWritten) {
                res.writeHead(404);
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

function mapQuery(app: express.Express, queryName: string, queryDefinition: QueryDefinition) {
    app.get(makePath(':entryId/' + queryName), async (req, res) => {
        ConsoleLogger.log('Server', `Query '${req.params.entryId}/${queryName}'...`);

        const request = createRequest('pdb', req.params.entryId, queryName, req.query);
        const writer = wrapResponse(request.responseFormat.isBinary ? 'result.bcif' : 'result.cif', res);
        writer.writeHeader(request.responseFormat.isBinary);
        await resolveRequest(request, writer);
        writer.end();
    });
}

export function initWebApi(app: express.Express) {
    for (const q of QueryList) {
        mapQuery(app, q.name, q.definition);
    }
}