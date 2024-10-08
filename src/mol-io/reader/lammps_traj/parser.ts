/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext, chunkedSubtask } from '../../../mol-task';
import { Tokenizer, TokenBuilder } from '../common/text/tokenizer';
import { ReaderResult as Result } from '../result';
import { TokenColumnProvider as TokenColumn } from '../common/text/column/token';
import { Column } from '../../../mol-data/db';

export interface LammpsFrame {
    count: number,
    atomId: Column<number>,
    moleculeId: Column<number>,
    atomType: Column<number>,
    x: Column<number>,
    y: Column<number>,
    z: Column<number>,
}

export interface LammpTrajectoryFile {
    frames: LammpsFrame[],
    times: number[],
    timeOffset: number,
    deltaTime: number
}

const { readLine, skipWhitespace, eatValue, eatLine, markStart } = Tokenizer;

function State(tokenizer: Tokenizer, runtimeCtx: RuntimeContext) {
    return {
        tokenizer,
        runtimeCtx,
    };
}
type State = ReturnType<typeof State>

async function handleAtoms(state: State, count: number): Promise<LammpsFrame> {
    const { tokenizer } = state;

    const atomId = TokenBuilder.create(tokenizer.data, count * 2);
    const moleculeType = TokenBuilder.create(tokenizer.data, count * 2);
    const atomType = TokenBuilder.create(tokenizer.data, count * 2);
    const x = TokenBuilder.create(tokenizer.data, count * 2);
    const y = TokenBuilder.create(tokenizer.data, count * 2);
    const z = TokenBuilder.create(tokenizer.data, count * 2);

    const { position } = tokenizer;
    tokenizer.position = position;

    const n = 6;

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
                    case 0: TokenBuilder.addUnchecked(atomId, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 1: TokenBuilder.addUnchecked(moleculeType, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 2: TokenBuilder.addUnchecked(atomType, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 3: TokenBuilder.addUnchecked(x, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 4: TokenBuilder.addUnchecked(y, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 5: TokenBuilder.addUnchecked(z, tokenizer.tokenStart, tokenizer.tokenEnd); break;
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
        atomId: TokenColumn(atomId)(Column.Schema.int),
        moleculeId: TokenColumn(moleculeType)(Column.Schema.int),
        atomType: TokenColumn(atomType)(Column.Schema.int),
        x: TokenColumn(x)(Column.Schema.float),
        y: TokenColumn(y)(Column.Schema.float),
        z: TokenColumn(z)(Column.Schema.float),
    };
}

async function parseInternal(data: string, ctx: RuntimeContext): Promise<Result<LammpTrajectoryFile>> {
    const tokenizer = Tokenizer(data);
    const state = State(tokenizer, ctx);
    const f: LammpTrajectoryFile = {
        frames: [],
        times: [],
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
            const frame: LammpsFrame = await handleAtoms(state, numAtoms);
            frames?.push(frame);
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

export function parseLammpTrajectory(data: string) {
    return Task.create<Result<LammpTrajectoryFile>>('Parse Lammp Trajectory', async ctx => {
        return await parseInternal(data, ctx);
    });
}