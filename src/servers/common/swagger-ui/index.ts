/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as express from 'express'
import * as fs from 'fs'
import { getAbsoluteFSPath } from 'swagger-ui-dist'
import { ServeStaticOptions } from 'serve-static';
import { interpolate } from 'mol-util/string';

export function swaggerUiAssetsHandler(options?: ServeStaticOptions) {
    const opts = options || {}
    opts.index = false
    return express.static(getAbsoluteFSPath(), opts)
}

function createHTML(swaggerUrl: string, apiPrefix: string) {
    const htmlTemplate = fs.readFileSync(`${__dirname}/indexTemplate.html`).toString()
    return interpolate(htmlTemplate, { swaggerUrl, apiPrefix })
}

export function swaggerUiIndexHandler(swaggerUrl: string, apiPrefix: string): express.Handler {
    const html = createHTML(swaggerUrl, apiPrefix)
    return (req: express.Request, res: express.Response) => {
        res.writeHead(200, { 'Content-Type': 'text/html; charset=utf-8' });
        res.end(html);
    }
}