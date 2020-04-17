/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as express from 'express';

import * as Api from './api';
import * as Data from './query/data-model';
import * as Coords from './algebra/coordinate';
import { ConsoleLogger } from '../../../mol-util/console-logger';
import { State } from './state';
import { LimitsConfig, ServerConfig } from '../config';
import { interpolate } from '../../../mol-util/string';
import { getSchema, shortcutIconLink } from './web-schema';
import { swaggerUiIndexHandler, swaggerUiAssetsHandler } from '../../common/swagger-ui';

export default function init(app: express.Express) {
    app.locals.mapFile = getMapFileFn();
    function makePath(p: string) {
        return `${ServerConfig.apiPrefix}/${p}`;
    }

    // Header
    app.get(makePath(':source/:id/?$'), (req, res) => getHeader(req, res));
    // Box /:src/:id/box/:a1,:a2,:a3/:b1,:b2,:b3?text=0|1&space=cartesian|fractional
    app.get(makePath(':source/:id/box/:a1,:a2,:a3/:b1,:b2,:b3/?'), (req, res) => queryBox(req, res, getQueryParams(req, false)));
    // Cell /:src/:id/cell/?text=0|1&space=cartesian|fractional
    app.get(makePath(':source/:id/cell/?'), (req, res) => queryBox(req, res, getQueryParams(req, true)));

    app.get(makePath('openapi.json'), (req, res) => {
        res.writeHead(200, {
            'Content-Type': 'application/json; charset=utf-8',
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Headers': 'X-Requested-With'
        });
        res.end(JSON.stringify(getSchema()));
    });

    app.use(makePath(''), swaggerUiAssetsHandler());
    app.get(makePath(''), swaggerUiIndexHandler({
        openapiJsonUrl: makePath('openapi.json'),
        apiPrefix: ServerConfig.apiPrefix,
        title: 'VolumeServer API',
        shortcutIconLink
    }));
}

function getMapFileFn() {
    const map = new Function('type', 'id', 'interpolate', [
        'id = id.toLowerCase()',
        'switch (type.toLowerCase()) {',
        ...ServerConfig.idMap.map(mapping => {
            const [type, path] = mapping;
            return `    case '${type}': return interpolate('${path}', { id });`;
        }),
        '    default: return void 0;',
        '}'
    ].join('\n'));
    return (type: string, id: string) => map(type, id, interpolate);
}

function wrapResponse(fn: string, res: express.Response) {
    return {
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
            return res.write(Buffer.from(data.buffer));
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
}

function getSourceInfo(req: express.Request) {
    return {
        filename: req.app.locals.mapFile(req.params.source, req.params.id),
        id: `${req.params.source}/${req.params.id}`
    };
}

function validateSourceAndId(req: express.Request, res: express.Response) {
    if (!req.params.source || req.params.source.length > 32 || !req.params.id || req.params.id.length > 32) {
        res.writeHead(404);
        res.end();
        ConsoleLogger.error(`Query Box`, 'Invalid source and/or id');
        return true;
    }
    return false;
}

async function getHeader(req: express.Request, res: express.Response) {
    if (validateSourceAndId(req, res)) {
        return;
    }

    let headerWritten = false;

    try {
        const { filename, id } = getSourceInfo(req);
        const header = await Api.getExtendedHeaderJson(filename, id);
        if (!header) {
            res.writeHead(404);
            return;
        }
        res.writeHead(200, {
            'Content-Type': 'application/json; charset=utf-8',
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Headers': 'X-Requested-With'
        });
        headerWritten = true;
        res.write(header);
    } catch (e) {
        ConsoleLogger.error(`Header ${req.params.source}/${req.params.id}`, e);
        if (!headerWritten) {
            res.writeHead(404);
        }
    } finally {
        res.end();
    }
}

function getQueryParams(req: express.Request, isCell: boolean): Data.QueryParams {
    const a = [+req.params.a1, +req.params.a2, +req.params.a3];
    const b = [+req.params.b1, +req.params.b2, +req.params.b3];

    const detail = Math.min(Math.max(0, (+req.query.detail) | 0), LimitsConfig.maxOutputSizeInVoxelCountByPrecisionLevel.length - 1);
    const isCartesian = (req.query.space as string || '').toLowerCase() !== 'fractional';

    const box: Data.QueryParamsBox = isCell
        ? { kind: 'Cell' }
        : (isCartesian
            ? { kind: 'Cartesian', a: Coords.cartesian(a[0], a[1], a[2]), b: Coords.cartesian(b[0], b[1], b[2]) }
            : { kind: 'Fractional', a: Coords.fractional(a[0], a[1], a[2]), b: Coords.fractional(b[0], b[1], b[2]) });

    const asBinary = (req.query.encoding as string || '').toLowerCase() !== 'cif';
    const sourceFilename = req.app.locals.mapFile(req.params.source, req.params.id)!;

    return {
        sourceFilename,
        sourceId: `${req.params.source}/${req.params.id}`,
        asBinary,
        box,
        detail
    };
}

async function queryBox(req: express.Request, res: express.Response, params: Data.QueryParams) {
    if (validateSourceAndId(req, res)) {
        return;
    }

    const outputFilename = Api.getOutputFilename(req.params.source, req.params.id, params);
    const response = wrapResponse(outputFilename, res);

    try {
        if (!params.sourceFilename) {
            response.do404();
            return;
        }

        let ok = await Api.queryBox(params, () => response);
        if (!ok) {
            response.do404();
            return;
        }
    } catch (e) {
        ConsoleLogger.error(`Query Box ${JSON.stringify(req.params || {})} | ${JSON.stringify(req.query || {})}`, e);
        response.do404();
    } finally {
        response.end();
        queryDone();
    }
}

function queryDone() {
    if (State.shutdownOnZeroPending) {
        process.exit(0);
    }
}