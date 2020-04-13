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
import * as zlib from 'zlib';

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

    endEntry() {
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
        return this.res.write(Buffer.from(data.buffer, data.byteOffset, data.byteLength));
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

        this.stream.pipe(this.res, { end: true });
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
        return !!this.stream.write(Buffer.from(data.buffer, data.byteOffset, data.byteLength));
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
    private file: fs.WriteStream | undefined = void 0;
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
        this.file = fs.createWriteStream(this.fn);
        this.opened = true;
    }

    writeBinary(data: Uint8Array) {
        this.open();
        this.file?.write(Buffer.from(data.buffer, data.byteOffset, data.byteLength));
        return true;
    }

    writeString(data: string) {
        this.open();
        this.file?.write(data);
        return true;
    }

    end() {
        if (!this.opened || this.ended) return;
        this.file?.end();
        this.ended = true;
    }

    constructor(private fn: string) {

    }
}

export class TarballFileResultWriter implements ResultWriter {
    private file: fs.WriteStream | undefined = void 0;
    private ended = false;
    private opened = false;
    private stream = zlib.createGzip({ level: this.gzipLevel, memLevel: 9, chunkSize: 16 * 16384 });
    private entrySize = 0;

    beginEntry(name: string, size: number) {
        const header = encodeTarHeader({ name, size });
        this.entrySize = size;
        this.stream.write(header);
    }

    endEntry() {
        const size = this.entrySize & 511;
        if (size) this.stream.write(END_OF_TAR.slice(0, 512 - size));
    }

    open() {
        if (this.opened) return;
        makeDir(path.dirname(this.fn));
        this.file = fs.createWriteStream(this.fn);
        this.stream.pipe(this.file, { end: true });

        this.opened = true;
    }

    writeBinary(data: Uint8Array) {
        this.open();
        this.stream.write(Buffer.from(data.buffer, data.byteOffset, data.byteLength));
        return true;
    }

    writeString(data: string) {
        this.open();
        this.stream.write(data);
        return true;
    }

    end() {
        if (!this.opened || this.ended) return;
        this.stream.write(END_OF_TAR);
        this.stream.end();
        this.ended = true;
    }

    constructor(private fn: string, private gzipLevel: number = 6) {

    }
}