/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, Model } from 'mol-model/structure';
import { PerformanceMonitor } from 'mol-util/performance-monitor';
import { Cache } from './cache';
import Config from '../config';
import CIF, { CifFrame } from 'mol-io/reader/cif'
import * as util from 'util'
import * as fs from 'fs'
import * as zlib from 'zlib'
import { Job } from './jobs';
import { ConsoleLogger } from 'mol-util/console-logger';
import { ModelPropertiesProvider } from '../property-provider';
import { parse_mmCIF } from 'mol-model-parsers/structure/mmcif';
import { ModelFormat } from 'mol-model-parsers/structure/format';

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

export async function createStructureWrapperFromJob(job: Job, propertyProvider: ModelPropertiesProvider | undefined, allowCache = true): Promise<StructureWrapper> {
    if (allowCache && Config.cacheParams.useCache) {
        const ret = StructureCache.get(job.key);
        if (ret) return ret;
    }
    const ret = await readStructureWrapper(job.key, job.sourceId, job.entryId, propertyProvider);
    if (allowCache && Config.cacheParams.useCache) {
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
        let input = await readFileAsync(filename)
        if (isGz) input = await unzipAsync(input);
        const data = new Uint8Array(input.byteLength);
        for (let i = 0; i < input.byteLength; i++) data[i] = input[i];
        return data;
    } else {
        if (isGz) {
            const data = await unzipAsync(await readFileAsync(filename));
            return data.toString('utf8');
        }
        return readFileAsync(filename, 'utf8');
    }
}

async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data);
    const parsed = await comp.run();
    if (parsed.isError) throw parsed;
    return parsed.result;
}

export async function readStructureWrapper(key: string, sourceId: string | '_local_', entryId: string, propertyProvider: ModelPropertiesProvider | undefined) {
    const filename = sourceId === '_local_' ? entryId : Config.mapFile(sourceId, entryId);
    if (!filename) throw new Error(`Cound not map '${key}' to a valid filename.`);
    if (!fs.existsSync(filename)) throw new Error(`Could not find source file for '${key}'.`);

    perf.start('read');
    let data;
    try {
        data = await readFile(filename);
    } catch (e) {
        ConsoleLogger.error(key, '' + e);
        throw new Error(`Could not read the file for '${key}' from disk.`);
    }

    perf.end('read');
    perf.start('parse');
    const frame = (await parseCif(data)).blocks[0];
    perf.end('parse');
    perf.start('createModel');
    const models = await parse_mmCIF(ModelFormat.mmCIF(frame)).run();
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
        isBinary: /\.bcif/.test(filename),
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