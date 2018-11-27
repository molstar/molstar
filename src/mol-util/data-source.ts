/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * Adapted from LiteMol
 */

import { Task, RuntimeContext } from 'mol-task';
import { utf8Read } from 'mol-io/common/utf8';

export enum DataCompressionMethod {
    None,
    Gzip
}

export interface AjaxGetParams {
    url: string,
    type: 'string' | 'binary',
    title?: string,
    compression?: DataCompressionMethod
}

export function readStringFromFile(file: File) {
    return <Task<string>>readFromFileInternal(file, false);
}

export function readUint8ArrayFromFile(file: File) {
    return <Task<Uint8Array>>readFromFileInternal(file, true);
}

export function readFromFile(file: File, type: 'string' | 'binary') {
    return <Task<Uint8Array | string>>readFromFileInternal(file, type === 'binary');
}

export function ajaxGetString(url: string, title?: string) {
    return <Task<string>>ajaxGetInternal(title, url, false, false);
}

export function ajaxGetUint8Array(url: string, title?: string) {
    return <Task<Uint8Array>>ajaxGetInternal(title, url, true, false);
}

export function ajaxGet(params: AjaxGetParams) {
    return <Task<string | Uint8Array>>ajaxGetInternal(params.title, params.url, params.type === 'binary', params.compression === DataCompressionMethod.Gzip);
}

function decompress(buffer: Uint8Array): Uint8Array {
    // TODO
    throw 'nyi';
    // const gzip = new LiteMolZlib.Gunzip(new Uint8Array(buffer));
    // return gzip.decompress();
}

async function processFile(ctx: RuntimeContext, asUint8Array: boolean, compressed: boolean, e: any) {
    const data = (e.target as FileReader).result;

    if (compressed) {
        await ctx.update('Decompressing...');

        const decompressed = decompress(new Uint8Array(data as ArrayBuffer));
        if (asUint8Array) {
            return decompressed;
        } else {
            return utf8Read(decompressed, 0, decompressed.length);
        }
    } else {
        return asUint8Array ? new Uint8Array(data as ArrayBuffer) : data as string;
    }
}

function readData(ctx: RuntimeContext, action: string, data: XMLHttpRequest | FileReader, asUint8Array: boolean): Promise<any> {
    return new Promise<any>((resolve, reject) => {
        data.onerror = (e: any) => {
            const error = (<FileReader>e.target).error;
            reject(error ? error : 'Failed.');
        };

        data.onabort = () => reject(Task.Aborted(''));

        data.onprogress = (e: ProgressEvent) => {
            if (e.lengthComputable) {
                ctx.update({ message: action, isIndeterminate: false, current: e.loaded, max: e.total });
            } else {
                ctx.update({ message: `${action} ${(e.loaded / 1024 / 1024).toFixed(2)} MB`, isIndeterminate: true });
            }
        }
        data.onload = (e: any) => resolve(e);
    });
}

function readFromFileInternal(file: File, asUint8Array: boolean): Task<string | Uint8Array> {
    let reader: FileReader | undefined = void 0;
    return Task.create('Read File', async ctx => {
        try {
            reader = new FileReader();
            const isCompressed = /\.gz$/i.test(file.name);

            if (isCompressed || asUint8Array) reader.readAsArrayBuffer(file);
            else reader.readAsBinaryString(file);

            ctx.update({ message: 'Opening file...', canAbort: true });
            const e = await readData(ctx, 'Reading...', reader, asUint8Array);
            const result = processFile(ctx, asUint8Array, isCompressed, e);
            return result;
        } finally {
            reader = void 0;
        }
    }, () => {
        if (reader) reader.abort();
    });
}

class RequestPool {
    private static pool: XMLHttpRequest[] = [];
    private static poolSize = 15;

    static get() {
        if (this.pool.length) {
            return this.pool.pop()!;
        }
        return new XMLHttpRequest();
    }

    static emptyFunc() { }

    static deposit(req: XMLHttpRequest) {
        if (this.pool.length < this.poolSize) {
            req.onabort = RequestPool.emptyFunc;
            req.onerror = RequestPool.emptyFunc;
            req.onload = RequestPool.emptyFunc;
            req.onprogress = RequestPool.emptyFunc;
            this.pool.push(req);
        }
    }
}

async function processAjax(ctx: RuntimeContext, asUint8Array: boolean, decompressGzip: boolean, e: any) {
    const req = (e.target as XMLHttpRequest);
    if (req.status >= 200 && req.status < 400) {
        if (asUint8Array) {
            const buff = new Uint8Array(e.target.response);
            RequestPool.deposit(e.target);

            if (decompressGzip) {
                return decompress(buff);
            } else {
                return buff;
            }
        }
        else {
            const text = e.target.responseText;
            RequestPool.deposit(e.target);
            return text;
        }
    } else {
        const status = req.statusText;
        RequestPool.deposit(e.target);
        throw status;
    }
}

function ajaxGetInternal(title: string | undefined, url: string, asUint8Array: boolean, decompressGzip: boolean): Task<string | Uint8Array> {
    let xhttp: XMLHttpRequest | undefined = void 0;
    return Task.create(title ? title : 'Download', async ctx => {
        try {
            if (!asUint8Array && decompressGzip) {
                throw 'Decompress is only available when downloading binary data.';
            }

            xhttp = RequestPool.get();

            xhttp.open('get', url, true);
            xhttp.responseType = asUint8Array ? 'arraybuffer' : 'text';
            xhttp.send();

            ctx.update({ message: 'Waiting for server...', canAbort: true });
            const e = await readData(ctx, 'Downloading...', xhttp, asUint8Array);
            const result = await processAjax(ctx, asUint8Array, decompressGzip, e)
            return result;
        } finally {
            xhttp = void 0;
        }
    }, () => {
        if (xhttp) xhttp.abort();
    });
}