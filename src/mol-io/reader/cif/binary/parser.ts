/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from '../data-model'
import * as Encoding from './encoding'
import Field from './field'
import Result from '../../result'
import decodeMsgPack from '../../../utils/msgpack/decode'
import Computation from '../../../../mol-base/computation'

function checkVersions(min: number[], current: number[]) {
    for (let i = 0; i < 2; i++) {
        if (min[i] > current[i]) return false;
    }
    return true;
}

function Category(data: Encoding.EncodedCategory): Data.Category {
    const map = Object.create(null);
    for (const col of data.columns) map[col.name] = col;
    return {
        rowCount: data.rowCount,
        getField(name) {
            const col = map[name];
            return col ? Field(col) : Data.DefaultUndefinedField(data.rowCount);
        }
    }
}

export default function parse(data: Uint8Array) {
    return Computation.create<Result<Data.File>>(async ctx => {
        const minVersion = [0, 3];

        try {
            const unpacked = decodeMsgPack(data) as Encoding.EncodedFile;
            if (!checkVersions(minVersion, unpacked.version.match(/(\d)\.(\d)\.\d/)!.slice(1).map(v => +v))) {
                return Result.error<Data.File>(`Unsupported format version. Current ${unpacked.version}, required ${minVersion.join('.')}.`);
            }
            const file = Data.File(unpacked.dataBlocks.map(block => {
                const cats = Object.create(null);
                for (const cat of block.categories) cats[cat.name] = Category(cat);
                return Data.Block(cats, block.header);
            }));
            return Result.success(file);
        } catch (e) {
            return Result.error<Data.File>('' + e);
        }
    })
}