/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, Model } from '../../../mol-model/structure';
import { PerformanceMonitor } from '../../../mol-util/performance-monitor';
import { Cache } from './cache';
import { ModelServerConfig as Config, mapSourceAndIdToFilename, ModelServerFetchFormats } from '../config';
import { CIF, CifFrame, CifBlock } from '../../../mol-io/reader/cif';
import * as util from 'util';
import * as fs from 'fs';
import * as zlib from 'zlib';
import { JobEntry } from './jobs';
import { ConsoleLogger } from '../../../mol-util/console-logger';
import { ModelPropertiesProvider } from '../property-provider';
import { trajectoryFromMmCIF } from '../../../mol-model-formats/structure/mmcif';
import { fetchRetry } from '../utils/fetch-retry';

require('util.promisify').shim();

export enum StructureSourceType {
    File,
    Cache
}

export interface StructureInfo {
    sourceType: StructureSourceType;
    readTime: number;
    parseTime: number;
    createModelTime: number;
    attachPropsTime: number;

    sourceId: string,
    entryId: string
}

export interface StructureWrapper {
    info: StructureInfo,
    isBinary: boolean,
    key: string,
    approximateSize: number,
    models: ArrayLike<Model>,
    modelMap: Map<number, Model>,
    structureModelMap: Map<number, Structure>,
    propertyProvider: ModelPropertiesProvider | undefined,
    cifFrame: CifFrame,
    cache: object
}

export async function createStructureWrapperFromJobEntry(entry: JobEntry, propertyProvider: ModelPropertiesProvider | undefined, allowCache = true): Promise<StructureWrapper> {
    if (allowCache && Config.cacheMaxSizeInBytes > 0) {
        const ret = StructureCache.get(entry.key);
        if (ret) return ret;
    }
    const ret = await readStructureWrapper(entry.key, entry.sourceId, entry.entryId, entry.job.id, propertyProvider);
    if (allowCache && Config.cacheMaxSizeInBytes > 0) {
        StructureCache.add(ret);
    }
    return ret;
}

export const StructureCache = new Cache<StructureWrapper>(s => s.key, s => s.approximateSize);
const perf = new PerformanceMonitor();

const readFileAsync = util.promisify(fs.readFile);
const unzipAsync = util.promisify<zlib.InputType, Buffer>(zlib.unzip);

async function readFile(filename: string) {
    const isGz = /\.gz$/i.test(filename);
    if (filename.match(/\.bcif/)) {
        let input = await readFileAsync(filename);
        if (isGz) input = await unzipAsync(input);
        const data = new Uint8Array(input.byteLength);
        for (let i = 0; i < input.byteLength; i++) data[i] = input[i];
        return { data, isBinary: true };
    } else {
        if (isGz) {
            const data = await unzipAsync(await readFileAsync(filename));
            return { data: data.toString('utf8'), isBinary: false };
        }
        return { data: await readFileAsync(filename, 'utf8'), isBinary: false };
    }
}

async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data);
    const parsed = await comp.run();
    if (parsed.isError) throw parsed;
    return parsed.result;
}

export async function readDataAndFrame(filename: string, key?: string): Promise<{ data: string | Uint8Array, frame: CifBlock, isBinary: boolean }> {
    perf.start('read');
    let data, isBinary;
    try {
        const read = await readFile(filename);
        data = read.data;
        isBinary = read.isBinary;
    } catch (e) {
        ConsoleLogger.error(key || filename, '' + e);
        throw new Error(`Could not read the file for '${key || filename}' from disk.`);
    }

    perf.end('read');
    perf.start('parse');
    const frame = (await parseCif(data)).blocks[0];
    perf.end('parse');

    return { data, frame, isBinary };
}

async function fetchDataAndFrame(jobId: string, uri: string, format: ModelServerFetchFormats, key?: string): Promise<{ data: string | Uint8Array, frame: CifBlock, isBinary: boolean }> {
    perf.start('read');
    const isBinary = format.startsWith('bcif');
    let data;
    try {
        ConsoleLogger.logId(jobId, 'Fetch', `${uri}`);
        const response = await fetchRetry(uri, 500, 3, () => ConsoleLogger.logId(jobId, 'Fetch', `Retrying to fetch '${uri}'`));

        if (format.endsWith('.gz')) {
            const input = await unzipAsync(await response.arrayBuffer());

            if (isBinary) {
                data = new Uint8Array(input.byteLength);
                for (let i = 0; i < input.byteLength; i++) data[i] = input[i];
            } else {
                data = input.toString('utf8');
            }
        } else {
            data = isBinary ? new Uint8Array(await response.arrayBuffer()) : await response.text();
        }
    } catch (e) {
        ConsoleLogger.error(key || uri, '' + e);
        throw new Error(`Could not fetch the file for '${key || uri}'.`);
    }

    perf.end('read');
    perf.start('parse');
    const frame = (await parseCif(data)).blocks[0];
    perf.end('parse');

    return { data, frame, isBinary };
}

function readOrFetch(jobId: string, key: string, sourceId: string | '_local_', entryId: string) {
    const mapped = sourceId === '_local_' ? [entryId] as const : mapSourceAndIdToFilename(sourceId, entryId);
    if (!mapped) throw new Error(`Cound not map '${key}' for a resource.`);

    const uri = mapped[0].toLowerCase();
    if (uri.startsWith('http://') || uri.startsWith('https://') || uri.startsWith('ftp://')) {
        return fetchDataAndFrame(jobId, mapped[0], (mapped[1] || 'cif').toLowerCase() as any, key);
    }

    if (!fs.existsSync(mapped[0])) throw new Error(`Could not find source file for '${key}'.`);
    return readDataAndFrame(mapped[0], key);
}

export async function readStructureWrapper(key: string, sourceId: string | '_local_', entryId: string, jobId: string | undefined, propertyProvider: ModelPropertiesProvider | undefined) {
    const { data, frame, isBinary } = await readOrFetch(jobId || '', key, sourceId, entryId);
    perf.start('createModel');
    const models = await trajectoryFromMmCIF(frame).run();
    perf.end('createModel');

    const modelMap = new Map<number, Model>();
    for (const m of models) {
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
            entryId
        },
        isBinary,
        key,
        approximateSize: typeof data === 'string' ? 2 * data.length : data.length,
        models,
        modelMap,
        structureModelMap: new Map(),
        cifFrame: frame,
        propertyProvider,
        cache: Object.create(null)
    };

    return ret;
}

export async function resolveStructure(wrapper: StructureWrapper, modelNum?: number) {
    if (typeof modelNum === 'undefined') modelNum = wrapper.models[0].modelNum;
    if (wrapper.structureModelMap.has(modelNum)) return wrapper.structureModelMap.get(modelNum)!;
    if (!wrapper.modelMap.has(modelNum)) {
        return void 0;
    }

    const model = wrapper.modelMap.get(modelNum)!;
    const structure = Structure.ofModel(model);
    if (wrapper.propertyProvider) {
        const modelProps = wrapper.propertyProvider(model, wrapper.cache);
        for (const p of modelProps) {
            await tryAttach(wrapper.key, p);
        }
    }
    return structure;
}

export async function resolveStructures(wrapper: StructureWrapper, modelNums?: number[]) {
    const ret: Structure[] = [];
    for (const n of modelNums || (wrapper.models as Model[]).map(m => m.modelNum)) {
        const s = await resolveStructure(wrapper, n);
        if (s) ret.push(s);
    }
    return ret;
}

async function tryAttach(key: string, promise: Promise<any>) {
    try {
        await promise;
    } catch (e) {
        ConsoleLogger.errorId(key, 'Custom prop:' + e);
    }
}