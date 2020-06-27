/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PdbFile } from './schema';
import { Task } from '../../../mol-task';
import { ReaderResult } from '../result';
import { Tokenizer } from '../common/text/tokenizer';

export function parsePDB(data: string, id?: string, isPdbqt = false): Task<ReaderResult<PdbFile>> {
    return Task.create('Parse PDB', async ctx => ReaderResult.success({
        lines: await Tokenizer.readAllLinesAsync(data, ctx),
        id,
        isPdbqt
    }));
}