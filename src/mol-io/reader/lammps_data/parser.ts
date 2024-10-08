/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Task, RuntimeContext, chunkedSubtask } from '../../../mol-task';
import { Tokenizer, TokenBuilder } from '../common/text/tokenizer';
import { ReaderResult as Result } from '../result';
import { TokenColumnProvider as TokenColumn } from '../common/text/column/token';
import { Column } from '../../../mol-data/db';

export interface LammpDataFile {
    readonly atoms: {
        readonly count: number
        readonly atomId: Column<number>
        readonly moleculeType: Column<number>
        readonly atomType: Column<number>
        readonly charge: Column<number>
        readonly x: Column<number>,
        readonly y: Column<number>,
        readonly z: Column<number>,
    }
    readonly bonds: {
        readonly count: number
        readonly bondId: Column<number>
        readonly bondType: Column<number>
        readonly atomIdA: Column<number>
        readonly atomIdB: Column<number>
    }
}

const { readLine, skipWhitespace, eatValue, eatLine, markStart } = Tokenizer;

const reWhitespace = /\s+/;

function State(tokenizer: Tokenizer, runtimeCtx: RuntimeContext) {
    return {
        tokenizer,
        runtimeCtx,
    };
}
type State = ReturnType<typeof State>

async function handleAtoms(state: State, count: number): Promise<LammpDataFile['atoms']> {
    const { tokenizer } = state;

    const atomId = TokenBuilder.create(tokenizer.data, count * 2);
    const moleculeType = TokenBuilder.create(tokenizer.data, count * 2);
    const atomType = TokenBuilder.create(tokenizer.data, count * 2);
    const charge = TokenBuilder.create(tokenizer.data, count * 2);
    const x = TokenBuilder.create(tokenizer.data, count * 2);
    const y = TokenBuilder.create(tokenizer.data, count * 2);
    const z = TokenBuilder.create(tokenizer.data, count * 2);

    const { position } = tokenizer;
    readLine(tokenizer).trim();
    tokenizer.position = position;

    const n = 7;

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
                    case 3: TokenBuilder.addUnchecked(charge, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 4: TokenBuilder.addUnchecked(x, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 5: TokenBuilder.addUnchecked(y, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 6: TokenBuilder.addUnchecked(z, tokenizer.tokenStart, tokenizer.tokenEnd); break;
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
        moleculeType: TokenColumn(moleculeType)(Column.Schema.int),
        atomType: TokenColumn(atomType)(Column.Schema.int),
        charge: TokenColumn(charge)(Column.Schema.float),
        x: TokenColumn(x)(Column.Schema.float),
        y: TokenColumn(y)(Column.Schema.float),
        z: TokenColumn(z)(Column.Schema.float),
    };
}

async function handleBonds(state: State, count: number): Promise<LammpDataFile['bonds']> {
    const { tokenizer } = state;

    const bondId = TokenBuilder.create(tokenizer.data, count * 2);
    const bondType = TokenBuilder.create(tokenizer.data, count * 2);
    const atomIdA = TokenBuilder.create(tokenizer.data, count * 2);
    const atomIdB = TokenBuilder.create(tokenizer.data, count * 2);

    const { length } = tokenizer;
    let bondsAlreadyRead = 0;
    await chunkedSubtask(state.runtimeCtx, 10, void 0, chunkSize => {
        const bondsToRead = Math.min(count - bondsAlreadyRead, chunkSize);
        for (let i = 0; i < bondsToRead; ++i) {
            for (let j = 0; j < 4; ++j) {
                skipWhitespace(tokenizer);
                markStart(tokenizer);
                eatValue(tokenizer);
                switch (j) {
                    case 0: TokenBuilder.addUnchecked(bondId, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 1: TokenBuilder.addUnchecked(bondType, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 2: TokenBuilder.addUnchecked(atomIdA, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 3: TokenBuilder.addUnchecked(atomIdB, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                }
            }
        }
        bondsAlreadyRead += bondsToRead;
        return bondsToRead;
    }, ctx => ctx.update({ message: 'Parsing...', current: tokenizer.position, max: length }));

    return {
        count,
        bondId: TokenColumn(atomIdA)(Column.Schema.int),
        bondType: TokenColumn(atomIdB)(Column.Schema.int),
        atomIdA: TokenColumn(atomIdA)(Column.Schema.int),
        atomIdB: TokenColumn(atomIdB)(Column.Schema.int),
    };
}

async function parseInternal(data: string, ctx: RuntimeContext): Promise<Result<LammpDataFile>> {
    const tokenizer = Tokenizer(data);
    const state = State(tokenizer, ctx);

    let atoms = undefined as LammpDataFile['atoms'] | undefined;
    let bonds = undefined as LammpDataFile['bonds'] | undefined;
    let numAtoms = 0;
    let numBonds = 0;
    while (tokenizer.tokenEnd < tokenizer.length) {
        const line = readLine(state.tokenizer).trim();
        if (line.includes('atoms')) {
            numAtoms = parseInt(line.split(reWhitespace)[0]);
        } else if (line.includes('bonds')) {
            numBonds = parseInt(line.split(reWhitespace)[0]);
        } else if (line.includes('Masses')) {
            // const numAtoms = parseInt(line.split(reWhitespace)[0]);
            // atoms = await handleAtoms(state, numAtoms);
        } else if (line.includes('Atoms')) {
            atoms = await handleAtoms(state, numAtoms);
        } else if (line.includes('Bonds')) {
            bonds = await handleBonds(state, numBonds);
        }
    }

    if (atoms === undefined) {
        return Result.error('no atoms data');
    }

    if (bonds === undefined) {
        return Result.error('no bonds data');
    }

    const result: LammpDataFile = {
        atoms,
        bonds
    };
    return Result.success(result);
}

export function parseLammpData(data: string) {
    return Task.create<Result<LammpDataFile>>('Parse LammpData', async ctx => {
        return await parseInternal(data, ctx);
    });
}