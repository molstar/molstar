/**
 * Copyright (c) 2021-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Jason Pattle <jpattle@exscientia.co.uk>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 */

import { Column } from '../../../mol-data/db';
import { MolFile } from '../mol/parser';
import { Tokenizer, TokenBuilder, Tokens } from '../common/text/tokenizer';
import { TokenColumnProvider as TokenColumn } from '../common/text/column/token';

export function isV3(
    versionLine: string
): boolean {
    return versionLine.trim().endsWith('V3000');
}

export function handleCountsV3(
    tokenizer: Tokenizer
): { atomCount: number, bondCount: number } {
    const atomCount = TokenBuilder.create(tokenizer.data, 1);
    const bondCount = TokenBuilder.create(tokenizer.data, 1);

    Tokenizer.eatLine(tokenizer); // BEGIN CTAB
    skipSingleValue(tokenizer); // M
    skipSingleValue(tokenizer); // V30
    skipSingleValue(tokenizer); // COUNTS

    addSingleValue(tokenizer, atomCount);
    addSingleValue(tokenizer, bondCount);
    Tokenizer.eatLine(tokenizer);

    return {
        atomCount: TokenColumn(atomCount)(Column.Schema.int).value(0),
        bondCount: TokenColumn(bondCount)(Column.Schema.int).value(0)
    };
}

export function handleAtomsV3(
    tokenizer: Tokenizer,
    atomCount: number
): MolFile['atoms'] {
    const x = TokenBuilder.create(tokenizer.data, atomCount * 2);
    const y = TokenBuilder.create(tokenizer.data, atomCount * 2);
    const z = TokenBuilder.create(tokenizer.data, atomCount * 2);
    const type_symbol = TokenBuilder.create(tokenizer.data, atomCount * 2);

    for (let i = 0; i < atomCount; ++i) {
        Tokenizer.markLine(tokenizer);
        skipSingleValue(tokenizer); // M
        skipSingleValue(tokenizer); // V30
        skipSingleValue(tokenizer); // Index

        const { position } = tokenizer;
        addSingleValue(tokenizer, type_symbol);
        addSingleValue(tokenizer, x);
        addSingleValue(tokenizer, y);
        addSingleValue(tokenizer, z);
        tokenizer.position = position;
    }
    Tokenizer.eatLine(tokenizer); // Previous Line
    Tokenizer.eatLine(tokenizer); // END ATOM

    return {
        count: atomCount,
        x: TokenColumn(x)(Column.Schema.float),
        y: TokenColumn(y)(Column.Schema.float),
        z: TokenColumn(z)(Column.Schema.float),
        type_symbol: TokenColumn(type_symbol)(Column.Schema.str),
        /* No support for formal charge parsing in V3000 molfiles at the moment,
        so all charges default to 0.*/
        formal_charge: Column.ofConst(0, atomCount, Column.Schema.int)
    };
}

export function handleBondsV3(
    tokenizer: Tokenizer,
    bondCount: number
): MolFile['bonds'] {
    const atomIdxA = TokenBuilder.create(tokenizer.data, bondCount * 2);
    const atomIdxB = TokenBuilder.create(tokenizer.data, bondCount * 2);
    const order = TokenBuilder.create(tokenizer.data, bondCount * 2);

    for (let i = 0; i < bondCount; ++i) {
        Tokenizer.markLine(tokenizer);
        skipSingleValue(tokenizer); // M
        skipSingleValue(tokenizer); // V30
        skipSingleValue(tokenizer); // Index

        const { position } = tokenizer;
        addSingleValue(tokenizer, order);
        addSingleValue(tokenizer, atomIdxA);
        addSingleValue(tokenizer, atomIdxB);
        tokenizer.position = position;
    }
    Tokenizer.eatLine(tokenizer); // Previous Line
    Tokenizer.eatLine(tokenizer); // END BOND

    return {
        count: bondCount,
        atomIdxA: TokenColumn(atomIdxA)(Column.Schema.float),
        atomIdxB: TokenColumn(atomIdxB)(Column.Schema.float),
        order: TokenColumn(order)(Column.Schema.float),
    };
}

function skipSingleValue(tokenizer: Tokenizer) {
    Tokenizer.skipWhitespace(tokenizer);
    Tokenizer.eatValue(tokenizer);
}

function addSingleValue(tokenizer: Tokenizer, tokens: Tokens) {
    const { position: valueStart } = tokenizer;
    Tokenizer.skipWhitespace(tokenizer);
    Tokenizer.eatValue(tokenizer);
    Tokenizer.trim(tokenizer, valueStart, tokenizer.position);
    TokenBuilder.addUnchecked(tokens, tokenizer.tokenStart, tokenizer.tokenEnd);
}
