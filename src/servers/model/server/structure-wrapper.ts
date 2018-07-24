/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, Model, Format } from 'mol-model/structure';
import { PerformanceMonitor } from 'mol-util/performance-monitor';
import { Cache } from './cache';
import Config from '../config';
import CIF from 'mol-io/reader/cif'
import * as util from 'util'
import * as fs from 'fs'
import * as zlib from 'zlib'
import { Job } from './jobs';
import { ConsoleLogger } from 'mol-util/console-logger';
import { attachModelProperties } from '../properties';

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

export class StructureWrapper {
    info: StructureInfo;

    key: string;
    approximateSize: number;
    structure: Structure;
}

export async function getStructure(job: Job): Promise<StructureWrapper> {
    if (Config.cacheParams.useCache) {
        const ret = StructureCache.get(job.key);
        if (ret) return ret;
    }
    const ret = await readStructure(job.key, job.sourceId, job.entryId);
    if (Config.cacheParams.useCache) {
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

async function readStructure(key: string, sourceId: string, entryId: string) {
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
    const models = await Model.create(Format.mmCIF(frame)).run();
    perf.end('createModel');

    perf.start('attachProps');
    await attachModelProperties(models[0]);
    perf.end('attachProps');

    const structure = Structure.ofModel(models[0]);

    const ret: StructureWrapper = {
        info: {
            sourceType: StructureSourceType.File,
            readTime: perf.time('read'),
            parseTime: perf.time('parse'),
            createModelTime: perf.time('createModel'),
            attachPropsTime: perf.time('attachProps'),
            sourceId,
            entryId
        },
        key,
        approximateSize: typeof data === 'string' ? 2 * data.length : data.length,
        structure
    };

    return ret;
}