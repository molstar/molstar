/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as express from 'express';
import * as fs from 'fs';
import * as path from 'path';
import { makeDir } from '../../../mol-util/make-dir';

export interface ResultWriter {
    beginEntry(name: string): void,
    endEntry(): void,

    writeBinary(data: Uint8Array): boolean,
    writeString(data: string): boolean,
    end(): void
}

export interface WebResutlWriter extends ResultWriter {
    doError(code?: number, message?: string): void,
    writeHeader(): void
}

export class SimpleResponseResultWriter implements WebResutlWriter {
    private ended = false;

    beginEntry(name: string) {
        throw new Error('Not supported');
    }

    endEntry() {
        throw new Error('Not supported');
    }

    doError(code = 404, message = 'Not Found.') {
        this.res.status(code).send(message);
        this.end();
    }

    writeHeader() {
        this.res.writeHead(200, {
            'Content-Type': this.isBinary ? 'application/octet-stream' : 'text/plain; charset=utf-8',
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Headers': 'X-Requested-With',
            'Content-Disposition': `inline; filename="${this.fn}"`
        });
    }

    writeBinary(data: Uint8Array) {
        return this.res.write(Buffer.from(data.buffer));
    }

    writeString(this: any, data: string) {
        return this.res.write(data);
    }

    end() {
        if (this.ended) return;
        this.res.end();
        this.ended = true;
    }

    constructor(private fn: string, private res: express.Response, private isBinary: boolean) {

    }
}

export class FileResultWriter implements ResultWriter {
    private file = 0;
    private ended = false;
    private opened = false;

    beginEntry(name: string) {
        throw new Error('Not supported');
    }

    endEntry() {
        throw new Error('Not supported');
    }

    open() {
        if (this.opened) return;
        makeDir(path.dirname(this.fn));
        this.file = fs.openSync(this.fn, 'w');
        this.opened = true;
    }

    writeBinary(data: Uint8Array) {
        this.open();
        fs.writeSync(this.file, Buffer.from(data.buffer));
        return true;
    }

    writeString(data: string) {
        this.open();
        fs.writeSync(this.file, data);
        return true;
    }

    end() {
        if (!this.opened || this.ended) return;
        fs.close(this.file, function () { });
        this.ended = true;
    }

    constructor(private fn: string) {

    }
}