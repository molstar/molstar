/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from LiteMol
 */

import { Task, RuntimeContext } from '../mol-task';

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

type DataType = 'json' | 'xml' | 'string' | 'binary'
type DataValue = 'string' | any | XMLDocument | Uint8Array
type DataResponse<T extends DataType> =
    T extends 'json' ? any :
        T extends 'xml' ? XMLDocument :
            T extends 'string' ? string :
                T extends 'binary' ? Uint8Array : never

export interface AjaxGetParams<T extends DataType = 'string'> {
    url: string,
    type?: T,
    title?: string,
    body?: string
}

export function readStringFromFile(file: File) {
    return readFromFileInternal(file, 'string');
}

export function readUint8ArrayFromFile(file: File) {
    return readFromFileInternal(file, 'binary');
}

export function readFromFile<T extends DataType>(file: File, type: T) {
    return readFromFileInternal(file, type);
}

export function ajaxGet(url: string): Task<DataValue>
export function ajaxGet<T extends DataType>(params: AjaxGetParams<T>): Task<DataResponse<T>>
export function ajaxGet<T extends DataType>(params: AjaxGetParams<T> | string) {
    if (typeof params === 'string') return ajaxGetInternal(params, params, 'string');
    return ajaxGetInternal(params.title, params.url, params.type || 'string', params.body);
}

export type AjaxTask = typeof ajaxGet

function isDone(data: XMLHttpRequest | FileReader) {
    if (data instanceof FileReader) {
        return data.readyState === FileReader.DONE
    } else if (data instanceof XMLHttpRequest) {
        return data.readyState === XMLHttpRequest.DONE
    }
    throw new Error('unknown data type')
}

function readData<T extends XMLHttpRequest | FileReader>(ctx: RuntimeContext, action: string, data: T): Promise<T> {
    return new Promise<T>((resolve, reject) => {
        // first check if data reading is already done
        if (isDone(data)) {
            const { error } = data as FileReader;
            if (error !== null) {
                reject(error ?? 'Failed.');
            } else {
                resolve(data);
            }
            return
        }

        data.onerror = (e: ProgressEvent) => {
            const { error } = e.target as FileReader;
            reject(error ?? 'Failed.');
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

function processFile<T extends DataType>(reader: FileReader, type: T): DataResponse<T> {
    const { result } = reader

    if (type === 'binary' && result instanceof ArrayBuffer) {
        return new Uint8Array(result) as DataResponse<T>
    } else if (type === 'string' && typeof result === 'string') {
        return result as DataResponse<T>
    } else if (type === 'xml' && typeof result === 'string') {
        const parser = new DOMParser();
        return parser.parseFromString(result, 'application/xml') as DataResponse<T>
    } else if (type === 'json' && typeof result === 'string') {
        return JSON.parse(result) as DataResponse<T>
    }
    throw new Error(`could not get requested response data '${type}'`)
}

function readFromFileInternal<T extends DataType>(file: File, type: T): Task<DataResponse<T>> {
    let reader: FileReader | undefined = void 0;
    return Task.create('Read File', async ctx => {
        try {
            reader = new FileReader();

            if (type === 'binary') reader.readAsArrayBuffer(file)
            else reader.readAsText(file)

            await ctx.update({ message: 'Opening file...', canAbort: true });
            const fileReader = await readData(ctx, 'Reading...', reader);

            await ctx.update({ message: 'Parsing file...', canAbort: false });
            return processFile(fileReader, type);
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

function processAjax<T extends DataType>(req: XMLHttpRequest, type: T): DataResponse<T> {
    if (req.status >= 200 && req.status < 400) {
        const { response } = req;
        RequestPool.deposit(req);

        if (type === 'binary' && response instanceof ArrayBuffer) {
            return new Uint8Array(response) as DataResponse<T>
        } else if (type === 'string' && typeof response === 'string') {
            return response as DataResponse<T>
        } else if (type === 'xml' && response instanceof XMLDocument) {
            return response as DataResponse<T>
        } else if (type === 'json' && typeof response === 'object') {
            return response as DataResponse<T>
        }
        throw new Error(`could not get requested response data '${type}'`)
    } else {
        const status = req.statusText;
        RequestPool.deposit(req);
        throw new Error(status);
    }
}

function getRequestResponseType(type: DataType): XMLHttpRequestResponseType {
    switch(type) {
        case 'json': return 'json'
        case 'xml': return 'document'
        case 'string': return 'text'
        case 'binary': return 'arraybuffer'
    }
}

function ajaxGetInternal<T extends DataType>(title: string | undefined, url: string, type: T, body?: string): Task<DataResponse<T>> {
    let xhttp: XMLHttpRequest | undefined = void 0;
    return Task.create(title ? title : 'Download', async ctx => {
        xhttp = RequestPool.get();

        xhttp.open(body ? 'post' : 'get', url, true);
        xhttp.responseType = getRequestResponseType(type);
        xhttp.send(body);

        await ctx.update({ message: 'Waiting for server...', canAbort: true });
        const req = await readData(ctx, 'Downloading...', xhttp);
        xhttp = void 0; // guard against reuse, help garbage collector

        await ctx.update({ message: 'Parsing response...', canAbort: false });
        const result = processAjax(req, type)

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