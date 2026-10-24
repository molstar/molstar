/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 *
 * Adapted from LiteMol
 */

import type { StringLike } from './string-like.js';
import { utf8Read, utf8ReadLong } from './utf8.js';
import { type RuntimeContext, Task } from '@molstar/core/task';
import { Asset, AssetManager } from './assets.js';
import { RUNNING_IN_NODEJS } from './nodejs-shims.js';
import { ungzip, unzip } from './zip/zip.js';

export enum DataCompressionMethod {
  None,
  Gzip,
  Zip,
}

export type DataType = 'json' | 'xml' | 'string' | 'binary' | 'zip';
export type DataValue = 'string' | any | XMLDocument | Uint8Array;
export type DataResponse<T extends DataType> = T extends 'json'
  ? any
  : T extends 'xml'
    ? XMLDocument
    : T extends 'string'
      ? StringLike
      : T extends 'binary'
        ? Uint8Array<ArrayBuffer>
        : T extends 'zip'
          ? { [k: string]: Uint8Array<ArrayBuffer> }
          : never;

export interface FetchDataParams<T extends DataType = 'string'> {
  url: string;
  type?: T;
  title?: string;
  headers?: Record<string, string>;
  body?: string;
  /**
   * Cancels the download, also while it waits for the server or on a stalled body, which the task's own abort (seen
   * between chunks) cannot; the task then rejects as aborted, with the signal's reason
   */
  signal?: AbortSignal;
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

/**
 * Downloads `url` (`fetch`; a `file://` URL in Node.js through `setFSModule`) and reads it as `type`, a GET or, with a
 * `body`, a POST.
 */
export function fetchData(url: string): Task<DataValue>;
export function fetchData<T extends DataType>(params: FetchDataParams<T>): Task<DataResponse<T>>;
export function fetchData<T extends DataType>(params: FetchDataParams<T> | string) {
  if (typeof params === 'string') return fetchDataInternal(params, params, 'string');
  return fetchDataInternal(
    params.title,
    params.url,
    params.type || 'string',
    params.body,
    params.headers,
    params.signal,
  );
}

export type FetchDataTask = typeof fetchData;

/** @deprecated Kept for compatibility: use `fetchData`, the same function (it has downloaded by `fetch` since v6). */
export const ajaxGet = fetchData;
/** @deprecated Kept for compatibility: use `FetchDataParams`. */
export type AjaxGetParams<T extends DataType = 'string'> = FetchDataParams<T>;
/** @deprecated Kept for compatibility: use `FetchDataTask`. */
export type AjaxTask = FetchDataTask;

function isDone(data: FileReader) {
  if (RUNNING_IN_NODEJS) throw new Error('`isDone` should not be used when running in Node.js'); // FileReader is not available in Node.js
  return data.readyState === FileReader.DONE;
}

function genericError(isDownload: boolean) {
  if (isDownload)
    return 'Failed to download data. Possible reasons: Resource is not available, or CORS is not allowed on the server.';
  return 'Failed to open file.';
}

function readData(ctx: RuntimeContext, action: string, data: FileReader): Promise<FileReader> {
  if (RUNNING_IN_NODEJS) throw new Error('`readData` should not be used when running in Node.js'); // FileReader is not available in Node.js
  return new Promise<FileReader>((resolve, reject) => {
    // first check if data reading is already done
    if (isDone(data)) {
      const { error } = data;
      if (error !== null && error !== undefined) {
        reject(error);
      } else {
        resolve(data);
      }
      return;
    }

    let hasError = false;

    data.onerror = (e: ProgressEvent) => {
      if (hasError) return;

      const { error } = e.target as FileReader;
      reject(error ?? genericError(false));
    };

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
    };

    data.onload = (e: ProgressEvent) => {
      resolve(data);
    };
  });
}

function getCompression(name: string) {
  return /\.gz$/i.test(name)
    ? DataCompressionMethod.Gzip
    : /\.zip$/i.test(name)
      ? DataCompressionMethod.Zip
      : DataCompressionMethod.None;
}

const reFilterPath = /^(__MACOSX|.DS_Store)/;

async function decompress(
  ctx: RuntimeContext,
  data: Uint8Array<ArrayBuffer>,
  compression: DataCompressionMethod,
): Promise<Uint8Array<ArrayBuffer>> {
  switch (compression) {
    case DataCompressionMethod.None:
      return data;
    case DataCompressionMethod.Gzip:
      return ungzip(ctx, data);
    case DataCompressionMethod.Zip:
      const parsed = await unzip(ctx, data.buffer);
      const names = Object.keys(parsed).filter((n) => !reFilterPath.test(n));
      if (names.length !== 1) throw new Error('can only decompress zip files with a single entry');
      return parsed[names[0]] as Uint8Array<ArrayBuffer>;
  }
}

async function processFile<T extends DataType>(
  ctx: RuntimeContext,
  fileContent: ArrayBuffer | null,
  type: T,
  compression: DataCompressionMethod,
): Promise<DataResponse<T>> {
  if (fileContent === null) throw new Error('no data given');
  return decodeData(ctx, new Uint8Array(fileContent), type, compression);
}

/**
 * Reads `data` as `type`, decompressing it first (`compression`, not for `'zip'`, which unzips). A `'binary'` read of
 * uncompressed data is `data` itself, not a copy.
 */
export async function decodeData<T extends DataType>(
  ctx: RuntimeContext,
  data: Uint8Array<ArrayBuffer>,
  type: T,
  compression: DataCompressionMethod = DataCompressionMethod.None,
): Promise<DataResponse<T>> {
  if (compression !== DataCompressionMethod.None && type !== 'zip') {
    // if type==='zip', data will be decompressed later
    data = await decompress(ctx, data, compression);
  }

  if (type === 'binary') {
    return data as DataResponse<T>;
  } else if (type === 'zip') {
    const buffer =
      data.byteOffset === 0 && data.byteLength === data.buffer.byteLength ? data.buffer : data.slice().buffer;
    return (await unzip(ctx, buffer)) as DataResponse<T>;
  } else if (type === 'string') {
    return utf8ReadLong(data) as DataResponse<T>;
  } else if (type === 'xml') {
    const parser = new DOMParser();
    return parser.parseFromString(utf8Read(data), 'application/xml') as DataResponse<T>;
  } else if (type === 'json') {
    return JSON.parse(utf8Read(data)) as DataResponse<T>;
  }
  throw new Error(`could not get requested response data '${type}'`);
}

function readFromFileInternal<T extends DataType>(file: File, type: T): Task<DataResponse<T>> {
  if (RUNNING_IN_NODEJS) {
    return readFromFileInternal_NodeJS(file, type);
  }
  let reader: FileReader | undefined = void 0;
  return Task.create(
    'Read File',
    async (ctx) => {
      try {
        await ctx.update({ message: 'Opening file...', canAbort: true });
        reader = new FileReader();
        reader.readAsArrayBuffer(file);
        const fileReader = await readData(ctx, 'Reading...', reader);
        const fileContent = fileReader.result as ArrayBuffer | null;

        await ctx.update({ message: 'Processing file...', canAbort: false });
        return await processFile(ctx, fileContent, type, getCompression(file.name));
      } finally {
        reader = void 0;
      }
    },
    () => {
      if (reader) reader.abort();
    },
  );
}
function readFromFileInternal_NodeJS<T extends DataType>(file: File, type: T): Task<DataResponse<T>> {
  return Task.create('Read File', async (ctx) => {
    await ctx.update({ message: 'Opening file...', canAbort: false });
    const fileContent = await file.arrayBuffer();

    await ctx.update({ message: 'Processing file...', canAbort: false });
    return await processFile(ctx, fileContent, type, getCompression(file.name));
  });
}

function fetchDataInternal<T extends DataType>(
  title: string | undefined,
  url: string,
  type: T,
  body?: string,
  headers?: Record<string, string>,
  signal?: AbortSignal,
): Task<DataResponse<T>> {
  // Node.js's `fetch` has no `file:` scheme
  if (RUNNING_IN_NODEJS && url.startsWith('file://')) return fetchDataInternal_file_NodeJS(title, url, type);

  return Task.create(title ? title : 'Download', async (ctx) => {
    // the task runtime sees an abort only in `ctx.update` and `ctx.shouldUpdate`, between chunks as with XHR; the
    // abort error leaves through the `finally` below, which cancels the request (no `onAbort`: it runs after `f` throws)
    const aborter = new AbortController();
    if (signal?.aborted) throw signalAborted(signal);
    // a listener, not `AbortSignal.any`: on Node.js 22 a long-lived caller signal keeps every combined signal alive
    const onCallerAbort = () => aborter.abort(signal!.reason);
    signal?.addEventListener('abort', onCallerAbort, { once: true });
    try {
      await ctx.update({ message: 'Waiting for server...', canAbort: true });
      let response: Response;
      try {
        response = await fetch(url, { signal: aborter.signal, method: body ? 'POST' : 'GET', body, headers });
      } catch (e) {
        throw downloadError(e, signal);
      }
      if (signal?.aborted) throw signalAborted(signal);
      if (!(response.status >= 200 && response.status < 400)) {
        throw new Error(`Download failed with status code ${response.status}`);
      }
      const data = await readResponseBody(ctx, 'Downloading...', response, signal);

      await ctx.update({ message: 'Parsing response...', canAbort: false });
      let result: DataResponse<T>;
      try {
        result = await processFile(ctx, data, type, DataCompressionMethod.None);
      } catch (e) {
        // XHR's `responseType: 'json'` gave null for a body that is not JSON, and callers fall back on it (`?? []`)
        if (!(type === 'json' && e instanceof SyntaxError)) throw e;
        result = null as DataResponse<T>;
      }
      if (signal?.aborted) throw signalAborted(signal);
      return result;
    } finally {
      signal?.removeEventListener('abort', onCallerAbort);
      // stops a body an error or an abort left unread; a no-op once it is read
      aborter.abort();
    }
  });
}

/**
 * The caller's signal aborts the task; otherwise a network failure. A browser hides why a request failed (XHR gave the
 * generic message); Node.js says why, as its path did before.
 */
function downloadError(e: unknown, signal: AbortSignal | undefined) {
  if (signal?.aborted) return signalAborted(signal);
  return RUNNING_IN_NODEJS ? e : new Error(genericError(true), { cause: e });
}

function signalAborted(signal: AbortSignal) {
  const { reason } = signal;
  return Task.Aborted(typeof reason === 'string' ? reason : reason instanceof Error ? reason.message : '');
}

/**
 * The bytes of a response body, with progress against `Content-Length`. A compressed body streams its decoded bytes,
 * which `Content-Length` does not count, so it reports bytes only, as does a body that outgrows the length: a
 * cross-origin response hides `Content-Encoding` but not `Content-Length`, so a compressed one shows as that.
 */
async function readResponseBody(
  ctx: RuntimeContext,
  action: string,
  response: Response,
  signal: AbortSignal | undefined,
): Promise<ArrayBuffer> {
  const encoding = response.headers.get('Content-Encoding');
  const length = encoding && encoding !== 'identity' ? -1 : contentLength(response);
  // a known length reads into one buffer, without a copy at the end
  let buffer = length >= 0 ? tryAllocate(length) : void 0;
  const chunks: Uint8Array[] = [];
  let loaded = 0;

  const reader = response.body?.getReader();
  while (reader) {
    let next: ReadableStreamReadResult<Uint8Array>;
    try {
      next = await reader.read();
    } catch (e) {
      throw downloadError(e, signal);
    }
    if (next.done) break;

    const chunk = next.value;
    if (buffer && loaded + chunk.length <= buffer.length) {
      buffer.set(chunk, loaded);
    } else {
      if (buffer) chunks.push(buffer.subarray(0, loaded));
      buffer = void 0;
      chunks.push(chunk);
    }
    loaded += chunk.length;

    if (ctx.shouldUpdate) {
      if (buffer) {
        await ctx.update({ message: action, isIndeterminate: false, current: loaded, max: length });
      } else {
        await ctx.update({ message: `${action} ${(loaded / 1024 / 1024).toFixed(2)} MB`, isIndeterminate: true });
      }
    }
  }

  if (buffer) return loaded === buffer.length ? buffer.buffer : buffer.slice(0, loaded).buffer;
  const out = new Uint8Array(loaded);
  for (let i = 0, offset = 0; i < chunks.length; offset += chunks[i].length, i++) out.set(chunks[i], offset);
  return out.buffer;
}

/** A buffer of `length` bytes, or none when that is more than can be allocated: a bogus length reads in chunks */
function tryAllocate(length: number) {
  try {
    return new Uint8Array(length);
  } catch (e) {
    if (e instanceof RangeError) return void 0;
    throw e;
  }
}

/** `Content-Length` in bytes, or -1 when the response has none */
function contentLength(response: Response) {
  const header = response.headers.get('Content-Length');
  const length = header === null ? NaN : Number(header);
  return Number.isSafeInteger(length) && length >= 0 ? length : -1;
}

// NOTE: a workaround for using this in Node.js
let _fs: typeof import('fs') | undefined = undefined;
function getFS() {
  if (!_fs) {
    throw new Error(
      "When running in Node.js and reading from files, call mol-util/data-source's setFSModule function first.",
    );
  }
  return _fs;
}

export function setFSModule(fs: typeof import('fs')) {
  _fs = fs;
}

function readFileAsync(filename: string): Promise<NonSharedBuffer> {
  return new Promise((resolve, reject) => {
    getFS().readFile(filename, (err, data) => {
      if (err) reject(err);
      else resolve(data);
    });
  });
}

/** Alternative implementation of fetchDataInternal for NodeJS for file:// protocol */
function fetchDataInternal_file_NodeJS<T extends DataType>(
  title: string | undefined,
  url: string,
  type: T,
): Task<DataResponse<T>> {
  if (!RUNNING_IN_NODEJS) throw new Error('This function should only be used when running in Node.js');
  if (!url.startsWith('file://')) throw new Error('This function is only for URLs with protocol file://');

  return Task.create(title ?? 'Download', async (ctx) => {
    const filename = url.substring('file://'.length);
    await ctx.update({ message: 'Loading file...', canAbort: false });
    const data = await readFileAsync(filename);

    await ctx.update({ message: 'Parsing response...', canAbort: false });
    const result = await processFile(ctx, data.buffer, type, DataCompressionMethod.None);
    return result;
  });
}

export type FetchDataManyEntry =
  | { kind: 'ok'; id: string; result: Asset.Wrapper<'string' | 'binary'> }
  | { kind: 'error'; id: string; error: any };
/** @deprecated Kept for compatibility: use `FetchDataManyEntry`. */
export type AjaxGetManyEntry = FetchDataManyEntry;

/**
 * A source of `fetchDataMany`. `shared`: a binary read of the held bytes (`AssetManager.ResolveOptions`), `decompress`:
 * gzip told by its bytes is inflated (`util/_spec/assets.md` §5).
 */
export interface FetchDataManySource {
  id: string;
  url: Asset.Url | string;
  isBinary?: boolean;
  canFail?: boolean;
  shared?: boolean;
  decompress?: boolean;
}

function resolveSource(ctx: RuntimeContext, assetManager: AssetManager, src: FetchDataManySource) {
  return AssetManager.resolveData(
    assetManager,
    Asset.getUrlAsset(assetManager, src.url),
    src.isBinary ? 'binary' : 'string',
    { shared: src.shared, decompress: src.decompress },
  ).runAsChild(ctx);
}

export async function fetchDataMany(
  ctx: RuntimeContext,
  assetManager: AssetManager,
  sources: FetchDataManySource[],
  maxConcurrency: number,
) {
  const len = sources.length;
  const slots: FetchDataManyEntry[] = new Array(sources.length);

  await ctx.update({ message: 'Downloading...', current: 0, max: len });
  let promises: Promise<FetchDataManyEntry & { index: number }>[] = [],
    promiseKeys: number[] = [];
  let currentSrc = 0;
  for (let _i = Math.min(len, maxConcurrency); currentSrc < _i; currentSrc++) {
    const current = sources[currentSrc];

    promises.push(wrapPromise(currentSrc, current.id, resolveSource(ctx, assetManager, current)));
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
      promises.push(wrapPromise(currentSrc, current.id, resolveSource(ctx, assetManager, current)));
      promiseKeys.push(currentSrc);
      currentSrc++;
    }
  }

  return slots;
}

function _filterRemoveIndex(this: number, _: any, i: number) {
  return this !== i;
}

async function wrapPromise(
  index: number,
  id: string,
  p: Promise<Asset.Wrapper<'string' | 'binary'>>,
): Promise<FetchDataManyEntry & { index: number }> {
  try {
    const result = await p;
    return { kind: 'ok', result, index, id };
  } catch (error) {
    return { kind: 'error', error, index, id };
  }
}

/** @deprecated Kept for compatibility: use `fetchDataMany`, the same function. */
export const ajaxGetMany = fetchDataMany;
