/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as express from 'express';
import * as fs from 'fs';
import * as path from 'path';
import { makeDir } from '../../../mol-util/make-dir';
import { encodeTarHeader, END_OF_TAR } from './tar';
import * as zlib from 'zlib'

export interface ResultWriter {
    beginEntry(name: string, size: number): void,
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
    private headerWritten = false;

    beginEntry(name: string) {
        throw new Error('Not supported');
    }

    async endEntry() {
        throw new Error('Not supported');
    }

    doError(code = 404, message = 'Not Found.') {
        if (!this.headerWritten) {
            this.headerWritten = true;
            this.res.status(code).send(message);
        }
        this.end();
    }

    writeHeader() {
        if (this.headerWritten) return;
        this.headerWritten = true;

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

export class TarballResponseResultWriter implements WebResutlWriter {
    private ended = false;
    private headerWritten = false;
    private stream = zlib.createGzip({ level: 6, memLevel: 9, chunkSize: 16 * 16384 });
    private entrySize = 0;


    beginEntry(name: string, size: number) {
        this.writeHeader();
        const header = encodeTarHeader({ name, size });
        this.entrySize = size;
        this.stream.write(header);
    }

    endEntry() {
        const size = this.entrySize & 511;
        if (size) this.stream.write(END_OF_TAR.slice(0, 512 - size));
    }

    doError(code = 404, message = 'Not Found.') {
        if (!this.headerWritten) {
            this.headerWritten = true;
            this.res.status(code).send(message);
        }
        this.end();
    }

    writeHeader() {
        if (this.headerWritten) return;

        this.stream.pipe(this.res);
        this.stream.on('end', () => this.res.end());

        this.headerWritten = true;
        this.res.writeHead(200, {
            'Content-Type': 'application/tar+gzip',
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Headers': 'X-Requested-With',
            'Content-Disposition': `inline; filename="${this.fn}"`
        });
    }

    writeBinary(data: Uint8Array) {
        this.writeHeader();
        return !!this.stream.write(Buffer.from(data.buffer));
    }

    writeString(data: string) {
        this.writeHeader();
        return !!this.stream.write(data);
    }

    end() {
        if (this.ended) return;
        this.ended = true;

        if (!this.headerWritten) {
            return;
        }

        this.stream.write(END_OF_TAR);
        this.stream.end();
    }

    constructor(private fn: string, private res: express.Response) {
    }
}

export class FileResultWriter implements ResultWriter {
    private file = 0;
    private ended = false;
    private opened = false;

    async beginEntry(name: string) {
        throw new Error('Not supported');
    }

    async endEntry() {
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