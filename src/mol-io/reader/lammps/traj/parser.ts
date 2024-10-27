/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Task, RuntimeContext, chunkedSubtask } from '../../../../mol-task';
import { Tokenizer, TokenBuilder } from '../../common/text/tokenizer';
import { ReaderResult as Result } from '../../result';
import { TokenColumnProvider as TokenColumn } from '../../common/text/column/token';
import { Column } from '../../../../mol-data/db';
import { LammpsFrame, LammpsTrajectoryFile } from '../schema';

const { readLine, skipWhitespace, eatValue, eatLine, markStart } = Tokenizer;

function State(tokenizer: Tokenizer, runtimeCtx: RuntimeContext) {
    return {
        tokenizer,
        runtimeCtx,
    };
}
type State = ReturnType<typeof State>

async function handleAtoms(state: State, count: number, parts: string[]): Promise<LammpsFrame> {
    const { tokenizer } = state;
    const columnIndexMap = Object.fromEntries(parts.map((colName, index) => [colName, index]));
    // declare column x, y, and z by check first caracter to 'x' or 'y' or 'z'
    // x,y,z = unscaled atom coordinates
    // xs,ys,zs = scaled atom coordinates this need the boundary box
    // xu,yu,zu = unwrapped atom coordinates
    // xsu,ysu,zsu = scaled unwrapped atom coordinates
    // ix,iy,iz = box image that the atom is in
    // how should we handle the different scenario ?
    const xCol = parts.findIndex(p => p[0] === 'x');
    const yCol = parts.findIndex(p => p[0] === 'y');
    const zCol = parts.findIndex(p => p[0] === 'z');
    // retrieve the atom type colum for x only
    const atomMode = parts[xCol]; // x,xs,xu,xsu
    const atomId = TokenBuilder.create(tokenizer.data, count * 2);
    const moleculeType = TokenBuilder.create(tokenizer.data, count * 2);
    const atomType = TokenBuilder.create(tokenizer.data, count * 2);
    const x = TokenBuilder.create(tokenizer.data, count * 2);
    const y = TokenBuilder.create(tokenizer.data, count * 2);
    const z = TokenBuilder.create(tokenizer.data, count * 2);

    const { position } = tokenizer;
    tokenizer.position = position;

    const n = parts.length;

    const { length } = tokenizer;
    let linesAlreadyRead = 0;
    await chunkedSubtask(state.runtimeCtx, 100000, void 0, chunkSize => {
        const linesToRead = Math.min(count - linesAlreadyRead, chunkSize);
        for (let i = 0; i < linesToRead; ++i) {
            for (let j = 0; j < n; ++j) {
                skipWhitespace(tokenizer);
                markStart(tokenizer);
                eatValue(tokenizer);
                switch (j) {
                    case columnIndexMap['id']: TokenBuilder.addUnchecked(atomId, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case columnIndexMap['mol']: TokenBuilder.addUnchecked(moleculeType, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case columnIndexMap['type']: TokenBuilder.addUnchecked(atomType, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case xCol: TokenBuilder.addUnchecked(x, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case yCol: TokenBuilder.addUnchecked(y, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case zCol: TokenBuilder.addUnchecked(z, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                }
            }
            // ignore any extra columns
            eatLine(tokenizer);
            markStart(tokenizer);
        }
        linesAlreadyRead += linesToRead;
        return linesToRead;
    }, ctx => ctx.update({ message: 'Parsing...', current: tokenizer.position, max: length }));

    return {
        count,
        atomMode: atomMode,
        atomId: TokenColumn(atomId)(Column.Schema.int),
        moleculeId: TokenColumn(moleculeType)(Column.Schema.int),
        atomType: TokenColumn(atomType)(Column.Schema.int),
        x: TokenColumn(x)(Column.Schema.float),
        y: TokenColumn(y)(Column.Schema.float),
        z: TokenColumn(z)(Column.Schema.float),
    };
}

/**
 * Possible Attributes from Lammps Dump
 * see https://docs.lammps.org/dump.html fro more details
 * possible attributes = id, mol, proc, procp1, type, element, mass,
 *                     x, y, z, xs, ys, zs, xu, yu, zu,
 *                     xsu, ysu, zsu, ix, iy, iz,
 *                     vx, vy, vz, fx, fy, fz,
 *                     q, mux, muy, muz, mu,
 *                     radius, diameter, omegax, omegay, omegaz,
 *                     angmomx, angmomy, angmomz, tqx, tqy, tqz,
 *                     c_ID, c_ID[I], f_ID, f_ID[I], v_name,
 *                     i_name, d_name, i2_name[I], d2_name[I]
 * ITEM: BOX BOUNDS xx yy zz
 * xlo xhi
 * ylo yhi
 * zlo zhi
 */
async function parseInternal(data: string, ctx: RuntimeContext): Promise<Result<LammpsTrajectoryFile>> {
    const tokenizer = Tokenizer(data);
    const state = State(tokenizer, ctx);
    const f: LammpsTrajectoryFile = {
        frames: [],
        times: [],
        bounds: [],
        timeOffset: 0.0,
        deltaTime: 0.0
    };
    const frames = f.frames;
    let numAtoms = 0;
    let timestep = 0;
    while (tokenizer.tokenEnd < tokenizer.length) {
        const line = readLine(state.tokenizer).trim();
        if (line.includes('ITEM: TIMESTEP')) {
            timestep = parseInt(readLine(state.tokenizer).trim());
            f.times.push(timestep);
        } else if (line.includes('ITEM: NUMBER OF ATOMS')) {
            numAtoms = parseInt(readLine(state.tokenizer).trim());
        } else if (line.includes('ITEM: ATOMS')) {
            // this line provide also the style of the output and will give the order of the columns
            const parts = line.split(' ').slice(2);
            const frame: LammpsFrame = await handleAtoms(state, numAtoms, parts);
            frames.push(frame);
        } else if (line.includes('ITEM: BOX BOUNDS')) {
            const tokens = line.split('ITEM: BOX BOUNDS ')[1].split(' ');
            // Periodicity of the box
            const px = tokens[0];
            const py = tokens[1];
            const pz = tokens[2];
            // the actual box bounds
            const xbound = readLine(state.tokenizer).trim().split(' ');
            const ybound = readLine(state.tokenizer).trim().split(' ');
            const zbound = readLine(state.tokenizer).trim().split(' ');
            const xlo = parseFloat(xbound[0]);
            const xhi = parseFloat(xbound[1]);
            const ylo = parseFloat(ybound[0]);
            const yhi = parseFloat(ybound[1]);
            const zlo = parseFloat(zbound[0]);
            const zhi = parseFloat(zbound[1]);
            f.bounds.push({
                lower: [xlo, ylo, zlo],
                length: [xhi - xlo, yhi - ylo, zhi - zlo],
                periodicity: [px, py, pz] });
        }
    }
    if (f.times.length >= 1) {
        f.timeOffset = f.times[0];
    }
    if (f.times.length >= 2) {
        f.deltaTime = f.times[1] - f.times[0];
    }
    return Result.success(f);
}

export function parseLammpsTrajectory(data: string) {
    return Task.create<Result<LammpsTrajectoryFile>>('Parse Lammp Trajectory', async ctx => {
        return await parseInternal(data, ctx);
    });
}