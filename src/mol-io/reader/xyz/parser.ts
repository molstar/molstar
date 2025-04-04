/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../../mol-data/db';
import { Task } from '../../../mol-task';
import { StringLike } from '../../common/string-like';
import { Tokenizer } from '../common/text/tokenizer';
import { ReaderResult as Result } from '../result';


export interface XyzFile {
    readonly molecules: {
        readonly comment: string,
        readonly count: number,
        readonly x: Column<number>,
        readonly y: Column<number>,
        readonly z: Column<number>,
        readonly type_symbol: Column<string>
    }[],
}

function handleMolecule(tokenizer: Tokenizer): XyzFile['molecules'][number] {
    let count = tokenizer.position >= tokenizer.data.length - 1 ? 0 : +Tokenizer.readLine(tokenizer);
    if (isNaN(count)) count = 0;

    const comment = Tokenizer.readLine(tokenizer);

    const x = new Float64Array(count);
    const y = new Float64Array(count);
    const z = new Float64Array(count);
    const type_symbol = new Array<string>(count);

    for (let i = 0; i < count; ++i) {
        const line = Tokenizer.readLineTrim(tokenizer);
        const fields = line.split(/\s+/g);
        type_symbol[i] = fields[0];
        x[i] = +fields[1];
        y[i] = +fields[2];
        z[i] = +fields[3];
    }

    return {
        count,
        comment,
        x: Column.ofFloatArray(x),
        y: Column.ofFloatArray(y),
        z: Column.ofFloatArray(z),
        type_symbol: Column.ofStringArray(type_symbol)
    };
}

function parseInternal(data: StringLike): Result<XyzFile> {
    const tokenizer = Tokenizer(data);

    const molecules: XyzFile['molecules'] = [];
    while (true) {
        const mol = handleMolecule(tokenizer);
        if (mol.count === 0) break;
        molecules.push(mol);
    }

    const result: XyzFile = { molecules };
    return Result.success(result);
}

export function parseXyz(data: StringLike) {
    return Task.create<Result<XyzFile>>('Parse Mol', async () => {
        return parseInternal(data);
    });
}