/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext, chunkedSubtask } from '../../../mol-task';
import { Tokenizer, TokenBuilder } from '../common/text/tokenizer';
import { ReaderResult as Result } from '../result';
import TokenColumn from '../common/text/column/token';
import { Column } from '../../../mol-data/db';

// http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html

export interface PsfFile {
    readonly id: string
    readonly title: string[]
    readonly atoms: {
        readonly count: number
        readonly atomId: Column<number>
        readonly segmentName: Column<string>
        readonly residueId: Column<number>
        readonly residueName: Column<string>
        readonly atomName: Column<string>
        readonly atomType: Column<string>
        readonly charge: Column<number>
        readonly mass: Column<number>
    }
    readonly bonds: {
        readonly count: number
        readonly atomIdA: Column<number>
        readonly atomIdB: Column<number>
    }
}

const { readLine, skipWhitespace, eatValue, eatLine, markStart } = Tokenizer;

const reWhitespace = /\s+/;
const reTitle = /(^\*|REMARK)*/;

function State(tokenizer: Tokenizer, runtimeCtx: RuntimeContext) {
    return {
        tokenizer,
        runtimeCtx,
    };
}
type State = ReturnType<typeof State>

async function handleAtoms(state: State, count: number): Promise<PsfFile['atoms']> {
    const { tokenizer } = state;

    const atomId = TokenBuilder.create(tokenizer.data, count * 2);
    const segmentName = TokenBuilder.create(tokenizer.data, count * 2);
    const residueId = TokenBuilder.create(tokenizer.data, count * 2);
    const residueName = TokenBuilder.create(tokenizer.data, count * 2);
    const atomName = TokenBuilder.create(tokenizer.data, count * 2);
    const atomType = TokenBuilder.create(tokenizer.data, count * 2);
    const charge = TokenBuilder.create(tokenizer.data, count * 2);
    const mass = TokenBuilder.create(tokenizer.data, count * 2);

    const { length } = tokenizer;
    let linesAlreadyRead = 0;
    await chunkedSubtask(state.runtimeCtx, 10, void 0, chunkSize => {
        const linesToRead = Math.min(count - linesAlreadyRead, chunkSize);
        for (let i = 0; i < linesToRead; ++i) {
            for (let j = 0; j < 8; ++j) {
                skipWhitespace(tokenizer);
                markStart(tokenizer);
                eatValue(tokenizer);
                switch (j) {
                    case 0: TokenBuilder.addUnchecked(atomId, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 1: TokenBuilder.addUnchecked(segmentName, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 2: TokenBuilder.addUnchecked(residueId, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 3: TokenBuilder.addUnchecked(residueName, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 4: TokenBuilder.addUnchecked(atomName, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 5: TokenBuilder.addUnchecked(atomType, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 6: TokenBuilder.addUnchecked(charge, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 7: TokenBuilder.addUnchecked(mass, tokenizer.tokenStart, tokenizer.tokenEnd); break;
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
        segmentName: TokenColumn(segmentName)(Column.Schema.str),
        residueId: TokenColumn(residueId)(Column.Schema.int),
        residueName: TokenColumn(residueName)(Column.Schema.str),
        atomName: TokenColumn(atomName)(Column.Schema.str),
        atomType: TokenColumn(atomType)(Column.Schema.str),
        charge: TokenColumn(charge)(Column.Schema.float),
        mass: TokenColumn(mass)(Column.Schema.float)
    };
}

async function handleBonds(state: State, count: number): Promise<PsfFile['bonds']> {
    const { tokenizer } = state;

    const atomIdA = TokenBuilder.create(tokenizer.data, count * 2);
    const atomIdB = TokenBuilder.create(tokenizer.data, count * 2);

    const { length } = tokenizer;
    let bondsAlreadyRead = 0;
    await chunkedSubtask(state.runtimeCtx, 10, void 0, chunkSize => {
        const bondsToRead = Math.min(count - bondsAlreadyRead, chunkSize);
        for (let i = 0; i < bondsToRead; ++i) {
            for (let j = 0; j < 2; ++j) {
                skipWhitespace(tokenizer);
                markStart(tokenizer);
                eatValue(tokenizer);
                switch (j) {
                    case 0: TokenBuilder.addUnchecked(atomIdA, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                    case 1: TokenBuilder.addUnchecked(atomIdB, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                }
            }
        }
        bondsAlreadyRead += bondsToRead;
        return bondsToRead;
    }, ctx => ctx.update({ message: 'Parsing...', current: tokenizer.position, max: length }));

    return {
        count,
        atomIdA: TokenColumn(atomIdA)(Column.Schema.int),
        atomIdB: TokenColumn(atomIdB)(Column.Schema.int),
    };
}

function parseTitle(state: State, count: number) {
    const title: string[] = [];
    for (let i = 0; i < count; ++i) {
        const line = readLine(state.tokenizer);
        title.push(line.replace(reTitle, '').trim());
    }
    return title;
}

async function parseInternal(data: string, ctx: RuntimeContext): Promise<Result<PsfFile>> {
    const tokenizer = Tokenizer(data);
    const state = State(tokenizer, ctx);

    let title = undefined as string[] | undefined;
    let atoms = undefined  as PsfFile['atoms'] | undefined;
    let bonds = undefined  as PsfFile['bonds'] | undefined;

    const id = readLine(state.tokenizer).trim();

    while(tokenizer.tokenEnd < tokenizer.length) {
        const line = readLine(state.tokenizer).trim();
        if (line.includes('!NTITLE')) {
            const numTitle = parseInt(line.split(reWhitespace)[0]);
            title = parseTitle(state, numTitle);
        } else if (line.includes('!NATOM')) {
            const numAtoms = parseInt(line.split(reWhitespace)[0]);
            atoms = await handleAtoms(state, numAtoms);
        } else if (line.includes('!NBOND')) {
            const numBonds = parseInt(line.split(reWhitespace)[0]);
            bonds = await handleBonds(state, numBonds);
            break; // TODO: don't break when the below are implemented
        } else if (line.includes('!NTHETA')) {
            // TODO
        } else if (line.includes('!NPHI')) {
            // TODO
        } else if (line.includes('!NIMPHI')) {
            // TODO
        } else if (line.includes('!NDON')) {
            // TODO
        } else if (line.includes('!NACC')) {
            // TODO
        } else if (line.includes('!NNB')) {
            // TODO
        } else if (line.includes('!NGRP NST2')) {
            // TODO
        } else if (line.includes('!MOLNT')) {
            // TODO
        } else if (line.includes('!NUMLP NUMLPH')) {
            // TODO
        } else if (line.includes('!NCRTERM')) {
            // TODO
        }
    }

    if (title === undefined) {
        title = [];
    }

    if (atoms === undefined) {
        return Result.error('no atoms data');
    }

    if (bonds === undefined) {
        return Result.error('no bonds data');
    }

    const result: PsfFile = {
        id,
        title,
        atoms,
        bonds
    };
    return Result.success(result);
}

export function parsePsf(data: string) {
    return Task.create<Result<PsfFile>>('Parse PSF', async ctx => {
        return await parseInternal(data, ctx);
    });
}