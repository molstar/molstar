/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as util from 'util';
import * as path from 'path';
import * as fs from 'fs';
import * as zlib from 'zlib';
import fetch from 'node-fetch';
require('util.promisify').shim();
const readFile = util.promisify(fs.readFile);
const writeFile = util.promisify(fs.writeFile);

import { Progress } from '../../mol-task';
import { Database } from '../../mol-data/db';
import { CIF } from '../../mol-io/reader/cif';
import { CifWriter } from '../../mol-io/writer/cif';
import { CCD_Schema } from '../../mol-io/reader/cif/schema/ccd';

export async function ensureAvailable(path: string, url: string, forceDownload = false) {
    if (forceDownload || !fs.existsSync(path)) {
        console.log(`downloading ${url}...`);
        const data = await fetch(url);
        if (!fs.existsSync(DATA_DIR)) {
            fs.mkdirSync(DATA_DIR);
        }
        if (url.endsWith('.gz')) {
            await writeFile(path, zlib.gunzipSync(await data.buffer()));
        } else {
            await writeFile(path, await data.text());
        }
        console.log(`done downloading ${url}`);
    }
}

export async function ensureDataAvailable(options: DataOptions) {
    await ensureAvailable(CCD_PATH, options.ccdUrl || CCD_URL, !!options.ccdUrl || options.forceDownload);
    await ensureAvailable(PVCD_PATH, options.pvcdUrl || PVCD_URL, !!options.pvcdUrl || options.forceDownload);
}

export async function readFileAsCollection<S extends Database.Schema>(path: string, schema: S) {
    const parsed = await parseCif(await readFile(path, 'utf8'));
    return CIF.toDatabaseCollection(schema, parsed.result);
}

export async function readCCD() {
    return readFileAsCollection(CCD_PATH, CCD_Schema);
}

export async function readPVCD() {
    return readFileAsCollection(PVCD_PATH, CCD_Schema);
}

async function parseCif(data: string | Uint8Array) {
    const comp = CIF.parse(data);
    console.time('parse cif');
    const parsed = await comp.run(p => console.log(Progress.format(p)), 250);
    console.timeEnd('parse cif');
    if (parsed.isError) throw parsed;
    return parsed;
}

export function getEncodedCif(name: string, database: Database<Database.Schema>, binary = false) {
    const encoder = CifWriter.createEncoder({ binary, encoderName: 'mol*' });
    CifWriter.Encoder.writeDatabase(encoder, name, database);
    return encoder.getData();
}

export type DataOptions = {
    ccdUrl?: string,
    pvcdUrl?: string,
    forceDownload?: boolean
}

export const DefaultDataOptions: DataOptions = {
    forceDownload: false
};

const DATA_DIR = path.join(__dirname, '..', '..', '..', '..', 'build/data');
const CCD_PATH = path.join(DATA_DIR, 'components.cif');
const PVCD_PATH = path.join(DATA_DIR, 'aa-variants-v1.cif');
const CCD_URL = 'https://files.wwpdb.org/pub/pdb/data/monomers/components.cif';
const PVCD_URL = 'https://files.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif';
