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
import { LammpsDataFile } from '../schema';
import { StringLike } from '../../../common/string-like';


const { readLine, skipWhitespace, eatValue, eatLine, markStart } = Tokenizer;

const reWhitespace = /\s+/;

function State(tokenizer: Tokenizer, runtimeCtx: RuntimeContext) {
    return {
        tokenizer,
        runtimeCtx,
    };
}
type State = ReturnType<typeof State>

async function handleAtoms(state: State, count: number, atom_style: 'full' | 'atomic' | 'bond'): Promise<LammpsDataFile['atoms']> {
    const { tokenizer } = state;
    // default atom style is atomic
    // depending on the atom style the number of columns can change
    const atomId = TokenBuilder.create(tokenizer.data, count * 2);
    const moleculeId = TokenBuilder.create(tokenizer.data, count * 2);
    const atomType = TokenBuilder.create(tokenizer.data, count * 2);
    const charge = TokenBuilder.create(tokenizer.data, count * 2);
    const x = TokenBuilder.create(tokenizer.data, count * 2);
    const y = TokenBuilder.create(tokenizer.data, count * 2);
    const z = TokenBuilder.create(tokenizer.data, count * 2);
    const columns = {
        full: [atomId, moleculeId, atomType, charge, x, y, z],
        atomic: [atomId, atomType, x, y, z],
        bond: [atomId, moleculeId, atomType, x, y, z],
    };
    const n = columns[atom_style].length;
    const { position } = tokenizer;
    readLine(tokenizer).trim();
    tokenizer.position = position;

    const { length } = tokenizer;
    let linesAlreadyRead = 0;

    await chunkedSubtask(state.runtimeCtx, 100000, void 0, chunkSize => {
        const linesToRead = Math.min(count - linesAlreadyRead, chunkSize);
        for (let i = 0; i < linesToRead; ++i) {
            for (let j = 0; j < n; ++j) {
                skipWhitespace(tokenizer);
                markStart(tokenizer);
                eatValue(tokenizer);
                const column = columns[atom_style][j];
                if (column) {
                    TokenBuilder.addUnchecked(column, tokenizer.tokenStart, tokenizer.tokenEnd);
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
        moleculeId: TokenColumn(moleculeId)(Column.Schema.int),
        atomType: TokenColumn(atomType)(Column.Schema.int),
        charge: TokenColumn(charge)(Column.Schema.float),
        x: TokenColumn(x)(Column.Schema.float),
        y: TokenColumn(y)(Column.Schema.float),
        z: TokenColumn(z)(Column.Schema.float),
    };
}

async function handleBonds(state: State, count: number): Promise<LammpsDataFile['bonds']> {
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
        bondId: TokenColumn(bondId)(Column.Schema.int),
        bondType: TokenColumn(bondType)(Column.Schema.int),
        atomIdA: TokenColumn(atomIdA)(Column.Schema.int),
        atomIdB: TokenColumn(atomIdB)(Column.Schema.int),
    };
}

const AtomStyles = ['full', 'atomic', 'bond'] as const;
type AtomStyle = typeof AtomStyles[number];

async function parseInternal(data: StringLike, ctx: RuntimeContext): Promise<Result<LammpsDataFile>> {
    const tokenizer = Tokenizer(data);
    const state = State(tokenizer, ctx);

    let atoms = undefined as LammpsDataFile['atoms'] | undefined;
    let bonds = undefined as LammpsDataFile['bonds'] | undefined;
    let numAtoms = 0;
    let numBonds = 0;
    let atom_style: AtomStyle = 'full';
    // full list of atom_style
    // https://docs.lammps.org/atom_style.html
    while (tokenizer.tokenEnd < tokenizer.length) {
        const line = readLine(state.tokenizer).trim();
        if (line.includes('atoms')) {
            numAtoms = parseInt(line.split(reWhitespace)[0]);
        } else if (line.includes('bonds')) {
            numBonds = parseInt(line.split(reWhitespace)[0]);
        } else if (line.includes('Masses')) {
            // TODO: support masses
        } else if (line.includes('Atoms')) {
            // usually atom style is indicated as a comment after Atoms. e.g. Atoms # full
            const parts = line.split('#');
            if (parts.length > 1) {
                const atomStyle = parts[1].trim();
                if (AtomStyles.includes(atomStyle as AtomStyle)) {
                    atom_style = atomStyle as AtomStyle;
                } else {
                    console.warn(`Unknown atom style: ${atomStyle}`);
                }
            }
            atoms = await handleAtoms(state, numAtoms, atom_style);
        } else if (line.includes('Bonds')) {
            bonds = await handleBonds(state, numBonds);
        }
    }

    if (atoms === undefined) {
        return Result.error('no atoms data');
    }

    if (bonds === undefined) {
        bonds = {
            count: 0,
            bondId: Column.ofIntArray([]),
            bondType: Column.ofIntArray([]),
            atomIdA: Column.ofIntArray([]),
            atomIdB: Column.ofIntArray([]),
        };
    }

    const result: LammpsDataFile = {
        atoms,
        bonds
    };
    return Result.success(result);
}

export function parseLammpsData(data: StringLike) {
    return Task.create<Result<LammpsDataFile>>('Parse LammpsData', async ctx => {
        return await parseInternal(data, ctx);
    });
}