/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from LiteMol
 */

import { Task, RuntimeContext } from '../mol-task';
import { utf8Read } from '../mol-io/common/utf8';
import { parseXml } from './xml-parser';
// polyfill XMLHttpRequest in node.js
const XHR = typeof document === 'undefined' ? require('xhr2') as {
    prototype: XMLHttpRequest;
    new(): XMLHttpRequest;
    readonly DONE: number;
    readonly HEADERS_RECEIVED: number;
    readonly LOADING: number;
    readonly OPENED: number;
    readonly UNSENT: number;
} : XMLHttpRequest

// export enum DataCompressionMethod {
//     None,
//     Gzip
// }

export interface AjaxGetParams<T extends 'string' | 'binary' | 'json' | 'xml' = 'string'> {
    url: string,
    type?: T,
    title?: string,
    // compression?: DataCompressionMethod
    body?: string
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

// TODO: support for no-referrer
export function ajaxGet(url: string): Task<string>
export function ajaxGet(params: AjaxGetParams<'string'>): Task<string>
export function ajaxGet(params: AjaxGetParams<'binary'>): Task<Uint8Array>
export function ajaxGet<T = any>(params: AjaxGetParams<'json' | 'xml'>): Task<T>
export function ajaxGet(params: AjaxGetParams<'string' | 'binary'>): Task<string | Uint8Array>
export function ajaxGet(params: AjaxGetParams<'string' | 'binary' | 'json' | 'xml'>): Task<string | Uint8Array | object>
export function ajaxGet(params: AjaxGetParams<'string' | 'binary' | 'json' | 'xml'> | string) {
    if (typeof params === 'string') return ajaxGetInternal(params, params, 'string', false);
    return ajaxGetInternal(params.title, params.url, params.type || 'string', false /* params.compression === DataCompressionMethod.Gzip */, params.body);
}

export type AjaxTask = typeof ajaxGet

function decompress(buffer: Uint8Array): Uint8Array {
    // TODO
    throw 'nyi';
    // const gzip = new LiteMolZlib.Gunzip(new Uint8Array(buffer));
    // return gzip.decompress();
}

async function processFile(ctx: RuntimeContext, asUint8Array: boolean, compressed: boolean, fileReader: FileReader) {
    const data = fileReader.result;

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

function isDone(data: XMLHttpRequest | FileReader) {
    if (data instanceof FileReader) {
        return data.readyState === FileReader.DONE
    } else if (data instanceof XMLHttpRequest) {
        return data.readyState === XMLHttpRequest.DONE
    }
    throw new Error('unknown data type')
}

function readData<T extends XMLHttpRequest | FileReader>(ctx: RuntimeContext, action: string, data: T, asUint8Array: boolean): Promise<T> {
    return new Promise<T>((resolve, reject) => {
        // first check if data reading is already done
        if (isDone(data)) {
            const error = (<FileReader>data).error;
            if (error) {
                reject((<FileReader>data).error || 'Failed.');
            } else {
                resolve(data);
            }
            return
        }

        data.onerror = (e: ProgressEvent) => {
            const error = (<FileReader>e.target).error;
            reject(error || 'Failed.');
        };

        let hasError = false;
        data.onprogress = (e: ProgressEvent) => {
            if (!ctx.shouldUpdate || hasError) return;

            try {
                if (e.lengthComputable) {
                    ctx.update({ message: action, isIndeterminate: false, current: e.loaded, max: e.total });
                } else {
                    ctx.update({ message: `${action} ${(e.loaded / 1024 / 1024).toFixed(2)} MB`, isIndeterminate: true });
                }
            } catch (e) {
                hasError = true;
                reject(e);
            }
        }

        data.onload = (e: ProgressEvent) => {
            resolve(data);
        }
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
            const fileReader = await readData(ctx, 'Reading...', reader, asUint8Array);
            const result = processFile(ctx, asUint8Array, isCompressed, fileReader);
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
        return new XHR();
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

async function processAjax(ctx: RuntimeContext, asUint8Array: boolean, decompressGzip: boolean, req: XMLHttpRequest) {
    if (req.status >= 200 && req.status < 400) {
        if (asUint8Array === true) {
            const buff = new Uint8Array(req.response);
            RequestPool.deposit(req);

            if (decompressGzip) {
                return decompress(buff);
            } else {
                return buff;
            }
        } else {
            const text = req.responseText;
            RequestPool.deposit(req);
            return text;
        }
    } else {
        const status = req.statusText;
        RequestPool.deposit(req);
        throw status;
    }
}

function ajaxGetInternal(title: string | undefined, url: string, type: 'json' | 'xml' | 'string' | 'binary', decompressGzip: boolean, body?: string): Task<string | Uint8Array> {
    let xhttp: XMLHttpRequest | undefined = void 0;
    return Task.create(title ? title : 'Download', async ctx => {
        const asUint8Array = type === 'binary';
        if (!asUint8Array && decompressGzip) {
            throw 'Decompress is only available when downloading binary data.';
        }

        xhttp = RequestPool.get();

        xhttp.open(body ? 'post' : 'get', url, true);
        xhttp.responseType = asUint8Array ? 'arraybuffer' : 'text';
        xhttp.send(body);

        await ctx.update({ message: 'Waiting for server...', canAbort: true });
        const req = await readData(ctx, 'Downloading...', xhttp, asUint8Array);
        xhttp = void 0; // guard against reuse, help garbage collector
        const result = await processAjax(ctx, asUint8Array, decompressGzip, req)

        if (type === 'json') {
            await ctx.update({ message: 'Parsing JSON...', canAbort: false });
            return JSON.parse(result as string);
        } else if (type === 'xml') {
            await ctx.update({ message: 'Parsing XML...', canAbort: false });
            return parseXml(result as string);
        }

        return result;
    }, () => {
        if (xhttp) {
            xhttp.abort();
            xhttp = void 0; // guard against reuse, help garbage collector
        }
    });
}

export type AjaxGetManyEntry<T> = { kind: 'ok', id: string, result: T } | { kind: 'error', id: string, error: any }
export async function ajaxGetMany(ctx: RuntimeContext, sources: { id: string, url: string, isBinary?: boolean, canFail?: boolean }[], maxConcurrency: number) {
    const len = sources.length;
    const slots: AjaxGetManyEntry<string | Uint8Array>[] = new Array(sources.length);

    await ctx.update({ message: 'Downloading...', current: 0, max: len });
    let promises: Promise<AjaxGetManyEntry<any> & { index: number }>[] = [], promiseKeys: number[] = [];
    let currentSrc = 0;
    for (let _i = Math.min(len, maxConcurrency); currentSrc < _i; currentSrc++) {
        const current = sources[currentSrc];
        promises.push(wrapPromise(currentSrc, current.id, ajaxGet({ url: current.url, type: current.isBinary ? 'binary' : 'string' }).runAsChild(ctx)));
        promiseKeys.push(currentSrc);
    }

    let done = 0;
    while (promises.length > 0) {
        const r = await Promise.race(promises);
        const src = sources[r.index];
        const idx = promiseKeys.indexOf(r.index);
        done++;
        if (r.kind === 'error' && !src.canFail) {
            // TODO: cancel other downloads
            throw new Error(`${src.url}: ${r.error}`);
        }
        if (ctx.shouldUpdate) {
            await ctx.update({ message: 'Downloading...', current: done, max: len });
        }
        slots[r.index] = r;
        promises = promises.filter(_filterRemoveIndex, idx);
        promiseKeys = promiseKeys.filter(_filterRemoveIndex, idx);
        if (currentSrc < len) {
            const current = sources[currentSrc];
            promises.push(wrapPromise(currentSrc, current.id, ajaxGet({ url: current.url, type: current.isBinary ? 'binary' : 'string' }).runAsChild(ctx)));
            promiseKeys.push(currentSrc);
            currentSrc++;
        }
    }

    return slots;
}

function _filterRemoveIndex(this: number, _: any, i: number) {
    return this !== i;
}

async function wrapPromise<T>(index: number, id: string, p: Promise<T>): Promise<AjaxGetManyEntry<T> & { index: number }> {
    try {
        const result = await p;
        return { kind: 'ok', result, index, id };
    } catch (error) {
        return { kind: 'error', error, index, id }
    }
}