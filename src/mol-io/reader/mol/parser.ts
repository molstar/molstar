/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../../mol-data/db';
import { Task } from '../../../mol-task';
import TokenColumn from '../common/text/column/token';
import { TokenBuilder, Tokenizer } from '../common/text/tokenizer';
import { ReaderResult as Result } from '../result';

/** Subset of the MolFile V2000 format */
export interface MolFile {
    readonly title: string,
    readonly program: string,
    readonly comment: string,
    readonly atoms: {
        readonly count: number,
        readonly x: Column<number>,
        readonly y: Column<number>,
        readonly z: Column<number>,
        readonly type_symbol: Column<string>
    },
    readonly bonds: {
        readonly count: number
        readonly atomIdxA: Column<number>,
        readonly atomIdxB: Column<number>,
        readonly order: Column<number>
    }
}

function handleAtoms(tokenizer: Tokenizer, count: number): MolFile['atoms'] {
    const x = TokenBuilder.create(tokenizer.data, count * 2);
    const y = TokenBuilder.create(tokenizer.data, count * 2);
    const z = TokenBuilder.create(tokenizer.data, count * 2);
    const type_symbol = TokenBuilder.create(tokenizer.data, count * 2);

    for (let i = 0; i < count; ++i) {
        Tokenizer.markLine(tokenizer);
        const { tokenStart: s, position } = tokenizer;
        Tokenizer.trim(tokenizer, s, s + 10);
        TokenBuilder.addUnchecked(x, tokenizer.tokenStart, tokenizer.tokenEnd);
        Tokenizer.trim(tokenizer, s + 10, s + 20);
        TokenBuilder.addUnchecked(y, tokenizer.tokenStart, tokenizer.tokenEnd);
        Tokenizer.trim(tokenizer, s + 20, s + 30);
        TokenBuilder.addUnchecked(z, tokenizer.tokenStart, tokenizer.tokenEnd);
        Tokenizer.trim(tokenizer, s + 31, s + 34);
        TokenBuilder.addUnchecked(type_symbol, tokenizer.tokenStart, tokenizer.tokenEnd);
        tokenizer.position = position;
    }

    return {
        count,
        x: TokenColumn(x)(Column.Schema.float),
        y: TokenColumn(y)(Column.Schema.float),
        z: TokenColumn(z)(Column.Schema.float),
        type_symbol: TokenColumn(type_symbol)(Column.Schema.str)
    };
}

function handleBonds(tokenizer: Tokenizer, count: number): MolFile['bonds'] {
    const atomIdxA = TokenBuilder.create(tokenizer.data, count * 2);
    const atomIdxB = TokenBuilder.create(tokenizer.data, count * 2);
    const order = TokenBuilder.create(tokenizer.data, count * 2);

    for (let i = 0; i < count; ++i) {
        Tokenizer.markLine(tokenizer);
        const { tokenStart: s, position } = tokenizer;
        Tokenizer.trim(tokenizer, s, s + 3);
        TokenBuilder.addUnchecked(atomIdxA, tokenizer.tokenStart, tokenizer.tokenEnd);
        Tokenizer.trim(tokenizer, s + 3, s + 6);
        TokenBuilder.addUnchecked(atomIdxB, tokenizer.tokenStart, tokenizer.tokenEnd);
        Tokenizer.trim(tokenizer, s + 6, s + 9);
        TokenBuilder.addUnchecked(order, tokenizer.tokenStart, tokenizer.tokenEnd);
        tokenizer.position = position;
    }

    return {
        count,
        atomIdxA: TokenColumn(atomIdxA)(Column.Schema.int),
        atomIdxB: TokenColumn(atomIdxB)(Column.Schema.int),
        order: TokenColumn(order)(Column.Schema.int)
    };
}

function parseInternal(data: string): Result<MolFile> {
    const tokenizer = Tokenizer(data);

    const title = Tokenizer.readLine(tokenizer).trim();
    const program = Tokenizer.readLine(tokenizer).trim();
    const comment = Tokenizer.readLine(tokenizer).trim();

    const counts = Tokenizer.readLine(tokenizer);

    const atomCount = +counts.substr(0, 3), bondCount = +counts.substr(3, 3);

    const atoms = handleAtoms(tokenizer, atomCount);
    const bonds = handleBonds(tokenizer, bondCount);

    const result: MolFile = {
        title,
        program,
        comment,
        atoms,
        bonds
    };
    return Result.success(result);
}

export function parseMol(data: string) {
    return Task.create<Result<MolFile>>('Parse Mol', async () => {
        return parseInternal(data);
    });
}