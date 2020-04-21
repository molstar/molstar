/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from NGL.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../../mol-math/linear-algebra';
import { chunkedSubtask, RuntimeContext, Task } from '../../../mol-task';
import { parseFloat as fastParseFloat } from '../common/text/number-parser';
import { Tokenizer } from '../common/text/tokenizer';
import { ReaderResult as Result } from '../result';

// http://apbs-pdb2pqr.readthedocs.io/en/latest/formats/opendx.html

export interface DxFile {
    header: DxFile.Header,
    values: Float64Array
}

export namespace DxFile {
    export interface Header {
        dim: Vec3,
        min: Vec3,
        h: Vec3
    }
}

function readHeader(tokenizer: Tokenizer): { header: DxFile.Header, headerByteCount: number } {
    const header: Partial<DxFile.Header> = { h: Vec3() };
    let headerByteCount = 0;
    let deltaLineCount = 0;

    const reWhitespace = /\s+/g;

    while (true) {
        const line = Tokenizer.readLine(tokenizer);
        let ls;

        if (line.startsWith('object 1')) {
            ls = line.split(reWhitespace);
            header.dim = Vec3.create(parseInt(ls[5]), parseInt(ls[6]), parseInt(ls[7]));
        } else if (line.startsWith('origin')) {
            ls = line.split(reWhitespace);
            header.min = Vec3.create(parseFloat(ls[1]), parseFloat(ls[2]), parseFloat(ls[3]));
        } else if (line.startsWith('delta')) {
            ls = line.split(reWhitespace);

            if (deltaLineCount === 0) {
                (header.h as any)[0] = parseFloat(ls[1]);
            } else if (deltaLineCount === 1) {
                (header.h as any)[1] = parseFloat(ls[2]);
            } else if (deltaLineCount === 2) {
                (header.h as any)[2] = parseFloat(ls[3]);
            }

            deltaLineCount += 1;
        } else if (line.startsWith('object 3')) {
            headerByteCount += line.length + 1;
            break;
        }

        headerByteCount += line.length + 1;
    }

    return { header: header as DxFile.Header, headerByteCount };
}

function readValuesText(ctx: RuntimeContext, tokenizer: Tokenizer, header: DxFile.Header) {
    const N = header.dim[0] * header.dim[1] * header.dim[2];
    const chunkSize = 100 * 100 * 100;
    const data = new Float64Array(N);
    let offset = 0;

    return chunkedSubtask(ctx, chunkSize, data, (count, data) => {
        const max = Math.min(N, offset + count);
        for (let i = offset; i < max; i++) {
            Tokenizer.skipWhitespace(tokenizer);
            tokenizer.tokenStart = tokenizer.position;
            Tokenizer.eatValue(tokenizer);
            data[i] = fastParseFloat(tokenizer.data, tokenizer.tokenStart, tokenizer.tokenEnd);
        }
        offset = max;
        return max === N ? 0 : chunkSize;
    }, (ctx, _, i) => ctx.update({ current: Math.min(i, N), max: N }));
}

export function parseDx(data: string | Uint8Array) {
    return Task.create<Result<DxFile>>('Parse Cube', async taskCtx => {
        await taskCtx.update('Reading header...');
        const tokenizer = Tokenizer(data as string);
        const { header } = readHeader(tokenizer);
        const values = await readValuesText(taskCtx, tokenizer, header);
        console.log(values);
        return Result.success({ header, values });
    });
}