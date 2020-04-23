/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from NGL.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../../mol-math/linear-algebra';
import { Tokenizer } from '../common/text/tokenizer';
import { Column } from '../../../mol-data/db';
import { Task, chunkedSubtask, RuntimeContext } from '../../../mol-task';
import { ReaderResult as Result } from '../result';
import { parseFloat as fastParseFloat } from '../common/text/number-parser';

// https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html

export interface CubeFile {
    name: string,
    header: CubeFile.Header,
    atoms: CubeFile.Atoms,
    values: Float64Array
}

export namespace CubeFile {
    export interface Header {
        orbitals: boolean,
        comment1: string,
        comment2: string,
        atomCount: number,
        origin: Vec3,
        dim: Vec3,
        basisX: Vec3,
        basisY: Vec3,
        basisZ: Vec3,
        dataSetIds: number[]
    }

    export interface Atoms {
        count: number,
        number: Column<number>,
        nuclearCharge: Column<number>,
        x: Column<number>,
        y: Column<number>,
        z: Column<number>
    }
}

const bohrToAngstromFactor = 0.529177210859;

function readHeader(tokenizer: Tokenizer) {
    const headerLines = Tokenizer.readLines(tokenizer, 6);
    const h = (k: number, l: number) => {
        const field = +headerLines[k].trim().split(/\s+/g)[ l ];
        return Number.isNaN(field) ? 0 : field;
    };
    const basis = (i: number) => {
        const n = h(i + 2, 0);
        const s = bohrToAngstromFactor;
        return [Math.abs(n), Vec3.create(h(i + 2, 1) * s, h(i + 2, 2) * s, h(i + 2, 3) * s), n] as const;
    };

    const comment1 = headerLines[0].trim();
    const comment2 = headerLines[1].trim();

    const [atomCount, origin, rawAtomCount] = basis(0);
    const [NVX, basisX] = basis(1);
    const [NVY, basisY] = basis(2);
    const [NVZ, basisZ] = basis(3);

    const atoms = readAtoms(tokenizer, atomCount, bohrToAngstromFactor);

    const dataSetIds: number[] = [];
    if (rawAtomCount >= 0) {
        let nVal = h(2, 4);
        if (nVal === 0) nVal = 1;
        for (let i = 0; i < nVal; i++) dataSetIds.push(i);
    } else {
        const counts = Tokenizer.readLine(tokenizer).trim().split(/\s+/g);
        for (let i = 0, _i = +counts[0]; i < _i; i++) dataSetIds.push(+counts[i + 1]);
    }

    const header: CubeFile.Header = { orbitals: rawAtomCount < 0, comment1, comment2, atomCount, origin, dim: Vec3.create(NVX, NVY, NVZ), basisX, basisY, basisZ, dataSetIds };
    return { header, atoms };
}

function readAtoms(tokenizer: Tokenizer, count: number, scaleFactor: number): CubeFile.Atoms {
    const number = new Int32Array(count);
    const value = new Float64Array(count);
    const x = new Float32Array(count);
    const y = new Float32Array(count);
    const z = new Float32Array(count);

    for (let i = 0; i < count; i++) {
        const fields = Tokenizer.readLine(tokenizer).trim().split(/\s+/g);
        number[i] = +fields[0];
        value[i] = +fields[1];
        x[i] = +fields[2] * scaleFactor;
        y[i] = +fields[3] * scaleFactor;
        z[i] = +fields[4] * scaleFactor;
    }

    return {
        count,
        number: Column.ofArray({ array: number, schema: Column.Schema.int }),
        nuclearCharge: Column.ofArray({ array: value, schema: Column.Schema.float }),
        x: Column.ofArray({ array: x, schema: Column.Schema.float }),
        y: Column.ofArray({ array: y, schema: Column.Schema.float }),
        z: Column.ofArray({ array: z, schema: Column.Schema.float })
    };
}

function readValues(ctx: RuntimeContext, tokenizer: Tokenizer, header: CubeFile.Header) {
    const N = header.dim[0] * header.dim[1] * header.dim[2] * header.dataSetIds.length;
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

export function parseCube(data: string, name: string) {
    return Task.create<Result<CubeFile>>('Parse Cube', async taskCtx => {
        await taskCtx.update('Reading header...');
        const tokenizer = Tokenizer(data);
        const { header, atoms } = readHeader(tokenizer);
        const values = await readValues(taskCtx, tokenizer, header);
        return Result.success({ header, atoms, values, name });
    });
}