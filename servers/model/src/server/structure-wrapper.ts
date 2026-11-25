/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { Structure, Model } from '@molstar/model/model/structure';
import { PerformanceMonitor } from '@molstar/core/util/performance-monitor';
import { Cache } from './cache.js';
import {
  ModelServerConfig as Config,
  mapSourceAndIdToFilename,
  ModelServerFetchFormats,
} from '@molstar/model-server/config';
import { CIF, type CifFrame, CifBlock } from '@molstar/io/reader/cif';
import * as util from 'util';
import * as fs from 'fs';
import * as zlib from 'zlib';
import { JobEntry } from './jobs.js';
import { ConsoleLogger } from '@molstar/core/util/console-logger';
import type { ModelPropertiesProvider } from '@molstar/model-server/property-provider';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { describeFetchError, fetchRetry } from '@molstar/model-server/utils/fetch-retry';
import { Task } from '@molstar/core/task';
import { HttpError } from './errors.js';

export enum StructureSourceType {
  File,
  Cache,
}

export interface StructureInfo {
  sourceType: StructureSourceType;
  readTime: number;
  parseTime: number;
  createModelTime: number;
  attachPropsTime: number;

  sourceId: string;
  entryId: string;
}

export interface StructureWrapper {
  info: StructureInfo;
  isBinary: boolean;
  key: string;
  /** Approximate memory the wrapper retains, in bytes (`estimateRetainedSize`). */
  approximateSize: number;
  models: ArrayLike<Model>;
  modelMap: Map<number, Model>;
  /** The structure of each model, built (and its properties attached) by the first query that needs it. */
  structureModelMap: Map<number, Promise<Structure>>;
  propertyProvider: ModelPropertiesProvider | undefined;
  cifFrame: CifFrame;
  cache: object;
}

/**
 * What a cached wrapper retains per byte of its (decompressed) source data, plus a fixed part per entry: the CIF
 * frame, the model and the model's structure after a query. Measured with `scripts/measure-retained.mjs` on RCSB
 * entries from 1cbs to 3j3q (bcif ×2.3–2.6, text ×2.7–2.8 on the large ones; small entries are dominated by the
 * fixed part); a `full` query, which decodes every column of a BinaryCIF, adds about 10% (3j3q: 249 → 272 MB).
 * `cacheMaxSizeInBytes` is a budget for these estimates.
 */
const RetainedBytesPerSourceByte = { binary: 2.5, text: 2.8 };
const RetainedBytesPerEntry = 2 * 1024 * 1024;

export function estimateRetainedSize(sourceLength: number, isBinary: boolean) {
  const perByte = isBinary ? RetainedBytesPerSourceByte.binary : RetainedBytesPerSourceByte.text;
  return Math.round(perByte * sourceLength) + RetainedBytesPerEntry;
}

/** Reads in progress by key: concurrent jobs for an entry that is not cached share one read. */
const pendingReads = new Map<string, Promise<StructureWrapper>>();

export async function createStructureWrapperFromJobEntry(
  entry: JobEntry,
  propertyProvider: ModelPropertiesProvider | undefined,
  allowCache = true,
): Promise<StructureWrapper> {
  const useCache = allowCache && Config.cacheMaxSizeInBytes > 0;
  if (useCache) {
    const ret = StructureCache.get(entry.key);
    if (ret) return ret;
  }

  let pending = pendingReads.get(entry.key);
  if (!pending) {
    pending = readStructureWrapper(entry.key, entry.sourceId, entry.entryId, entry.job.id, propertyProvider)
      .then((ret) => {
        if (useCache) StructureCache.add(ret);
        return ret;
      })
      .finally(() => pendingReads.delete(entry.key));
    pendingReads.set(entry.key, pending);
  }
  return pending;
}

export const StructureCache = new Cache<StructureWrapper>(
  (s) => s.key,
  (s) => s.approximateSize,
  () => ({ maxApproximateSizeInBytes: Config.cacheMaxSizeInBytes, entryTimeoutInMs: Config.cacheEntryTimeoutMs }),
);

const readFileAsync = fs.promises.readFile;
const unzipAsync = util.promisify<zlib.InputType, Buffer>(zlib.unzip);

/** A plain `Uint8Array` over the buffer's bytes, no copy: the bcif decoder copies what it keeps. */
function asUint8Array(buffer: Buffer | ArrayBuffer) {
  return buffer instanceof ArrayBuffer
    ? new Uint8Array(buffer)
    : new Uint8Array(buffer.buffer, buffer.byteOffset, buffer.byteLength);
}

async function readFile(filename: string) {
  const isGz = /\.gz$/i.test(filename);
  if (filename.match(/\.bcif/)) {
    let input = await readFileAsync(filename);
    if (isGz) input = (await unzipAsync(input)) as NonSharedBuffer;
    return { data: asUint8Array(input), isBinary: true };
  } else {
    if (isGz) {
      const data = await unzipAsync(await readFileAsync(filename));
      return { data: data.toString('utf8'), isBinary: false };
    }
    return { data: await readFileAsync(filename, 'utf8'), isBinary: false };
  }
}

async function parseCif(data: string | Uint8Array) {
  const comp = CIF.parse(data);
  const parsed = await comp.run();
  if (parsed.isError) throw parsed;
  return parsed.result;
}

export async function readDataAndFrame(
  filename: string,
  key?: string,
  perf = new PerformanceMonitor(),
): Promise<{ data: string | Uint8Array; frame: CifBlock; isBinary: boolean }> {
  perf.start('read');
  let data, isBinary;
  try {
    const read = await readFile(filename);
    data = read.data;
    isBinary = read.isBinary;
  } catch (e) {
    ConsoleLogger.error(key || filename, '' + e);
    throw new HttpError(404, `Could not read the file for '${key || filename}' from disk.`);
  }

  perf.end('read');
  perf.start('parse');
  const frame = (await parseCif(data)).blocks[0];
  perf.end('parse');

  return { data, frame, isBinary };
}

async function fetchDataAndFrame(
  jobId: string,
  uri: string,
  format: ModelServerFetchFormats,
  key: string | undefined,
  perf: PerformanceMonitor,
): Promise<{ data: string | Uint8Array; frame: CifBlock; isBinary: boolean }> {
  perf.start('read');
  const isBinary = format.startsWith('bcif');
  let data;
  try {
    ConsoleLogger.logId(jobId, 'Fetch', `${uri}`);
    const response = await fetchRetry(uri, Config.fetchTimeoutMs, 3, () =>
      ConsoleLogger.logId(jobId, 'Fetch', `Retrying to fetch '${uri}'`),
    );
    if (!response.ok) {
      void response.body?.cancel().catch(() => {});
      // source does not have the entry (4xx) or failed to serve it (5xx, after the retries)
      throw new HttpError(response.status < 500 ? 404 : 502, `HTTP ${response.status}`);
    }

    if (format.endsWith('.gz')) {
      const input = await unzipAsync(await response.arrayBuffer());
      data = isBinary ? asUint8Array(input) : input.toString('utf8');
    } else {
      data = isBinary ? asUint8Array(await response.arrayBuffer()) : await response.text();
    }
  } catch (e) {
    ConsoleLogger.error(key || uri, describeFetchError(e));
    throw new HttpError(e instanceof HttpError ? e.status : 502, `Could not fetch the file for '${key || uri}'.`);
  }

  perf.end('read');
  perf.start('parse');
  const frame = (await parseCif(data)).blocks[0];
  perf.end('parse');

  return { data, frame, isBinary };
}

function readOrFetch(
  jobId: string,
  key: string,
  sourceId: string | '_local_',
  entryId: string,
  perf: PerformanceMonitor,
) {
  const mapped = sourceId === '_local_' ? ([entryId] as const) : mapSourceAndIdToFilename(sourceId, entryId);
  if (!mapped) throw new HttpError(404, `Cound not map '${key}' for a resource.`);

  const uri = mapped[0].toLowerCase();
  if (uri.startsWith('http://') || uri.startsWith('https://') || uri.startsWith('ftp://') || uri.startsWith('gs://')) {
    return fetchDataAndFrame(jobId, mapped[0], (mapped[1] || 'cif').toLowerCase() as any, key, perf);
  }

  if (!fs.existsSync(mapped[0])) throw new HttpError(404, `Could not find source file for '${key}'.`);
  return readDataAndFrame(mapped[0], key, perf);
}

export async function readStructureWrapper(
  key: string,
  sourceId: string | '_local_',
  entryId: string,
  jobId: string | undefined,
  propertyProvider: ModelPropertiesProvider | undefined,
) {
  // one monitor per read: reads of concurrent jobs must not share timers
  const perf = new PerformanceMonitor();
  const { data, frame, isBinary } = await readOrFetch(jobId || '', key, sourceId, entryId, perf);
  perf.start('createModel');
  const trajectory = await trajectoryFromMmCIF(frame).run();
  perf.end('createModel');

  const models: Model[] = [];
  const modelMap = new Map<number, Model>();

  for (let i = 0; i < trajectory.frameCount; i++) {
    const m = await Task.resolveInContext(trajectory.getFrameAtIndex(i));
    models.push(m);
    modelMap.set(m.modelNum, m);
  }

  const ret: StructureWrapper = {
    info: {
      sourceType: StructureSourceType.File,
      readTime: perf.time('read'),
      parseTime: perf.time('parse'),
      createModelTime: perf.time('createModel'),
      attachPropsTime: 0, // perf.time('attachProps'),
      sourceId,
      entryId,
    },
    isBinary,
    key,
    approximateSize: estimateRetainedSize(data.length, isBinary),
    models,
    modelMap,
    structureModelMap: new Map(),
    cifFrame: frame,
    propertyProvider,
    cache: Object.create(null),
  };

  return ret;
}

export async function resolveStructure(wrapper: StructureWrapper, modelNum?: number) {
  if (typeof modelNum === 'undefined') modelNum = wrapper.models[0].modelNum;
  const model = wrapper.modelMap.get(modelNum);
  if (!model) return void 0;

  // promise is stored before the properties are attached, so concurrent queries share the work
  let structure = wrapper.structureModelMap.get(modelNum);
  if (!structure) {
    const created = createStructure(wrapper, model);
    wrapper.structureModelMap.set(modelNum, created);
    // failed builds aren't kept: the next query tries again
    created.catch(() => {
      if (wrapper.structureModelMap.get(modelNum) === created) wrapper.structureModelMap.delete(modelNum);
    });
    structure = created;
  }
  return structure;
}

async function createStructure(wrapper: StructureWrapper, model: Model) {
  const structure = Structure.ofModel(model);
  if (wrapper.propertyProvider) {
    let attaching: Promise<any>[] = [];
    try {
      attaching = wrapper.propertyProvider(model, wrapper.cache);
    } catch (e) {
      ConsoleLogger.errorId(wrapper.key, 'Custom prop:' + e);
    }
    // `allSettled`: the properties attach in parallel, and a failing one must not leave the others unhandled
    for (const r of await Promise.allSettled(attaching)) {
      if (r.status === 'rejected') ConsoleLogger.errorId(wrapper.key, 'Custom prop:' + r.reason);
    }
  }
  return structure;
}

export async function resolveStructures(wrapper: StructureWrapper, modelNums?: number[]) {
  const ret: Structure[] = [];
  for (const n of modelNums || (wrapper.models as Model[]).map((m) => m.modelNum)) {
    const s = await resolveStructure(wrapper, n);
    if (s) ret.push(s);
  }
  return ret;
}
