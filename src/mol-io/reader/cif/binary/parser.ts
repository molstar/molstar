/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from '../data-model';
import { EncodedCategory, EncodedFile } from '../../../common/binary-cif';
import Field from './field';
import { ReaderResult as Result } from '../../result';
import decodeMsgPack from '../../../common/msgpack/decode';
import { Task } from '../../../../mol-task';

function checkVersions(min: number[], current: number[]) {
    for (let i = 0; i < 2; i++) {
        if (min[i] > current[i]) return false;
    }
    return true;
}

function Category(data: EncodedCategory): Data.CifCategory {
    const map = Object.create(null);
    const cache = Object.create(null);
    for (const col of data.columns) map[col.name] = col;
    return {
        rowCount: data.rowCount,
        name: data.name.substr(1),
        fieldNames: data.columns.map(c => c.name),
        getField(name) {
            const col = map[name];
            if (!col) return void 0;
            if (!!cache[name]) return cache[name];
            cache[name] = Field(col);
            return cache[name];
        }
    };
}

export default function parse(data: Uint8Array) {
    return Task.create<Result<Data.CifFile>>('Parse BinaryCIF', async ctx => {
        const minVersion = [0, 3];

        try {
            const unpacked = decodeMsgPack(data) as EncodedFile;
            if (!checkVersions(minVersion, unpacked.version.match(/(\d)\.(\d)\.\d/)!.slice(1).map(v => +v))) {
                return Result.error<Data.CifFile>(`Unsupported format version. Current ${unpacked.version}, required ${minVersion.join('.')}.`);
            }
            const file = Data.CifFile(unpacked.dataBlocks.map(block => {
                const cats = Object.create(null);
                for (const cat of block.categories) cats[cat.name.substr(1)] = Category(cat);
                return Data.CifBlock(block.categories.map(c => c.name.substr(1)), cats, block.header);
            }));
            return Result.success(file);
        } catch (e) {
            return Result.error<Data.CifFile>('' + e);
        }
    });
}