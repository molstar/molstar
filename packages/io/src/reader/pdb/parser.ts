/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PdbFile } from './schema.js';
import { Task } from '@molstar/core/task';
import { ReaderResult } from '../result.js';
import { Tokenizer } from '../common/text/tokenizer.js';
import type { StringLike } from '@molstar/core/util/string-like';

export function parsePDB(data: StringLike, id?: string, variant?: 'pdb' | 'pdbqt' | 'pqr'): Task<ReaderResult<PdbFile>> {
    return Task.create('Parse PDB', async ctx => ReaderResult.success({
        lines: await Tokenizer.readAllLinesAsync(data, ctx),
        id,
        variant: variant || 'pdb',
    }));
}