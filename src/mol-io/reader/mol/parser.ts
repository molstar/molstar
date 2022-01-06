/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 */

import { Column } from '../../../mol-data/db';
import { Task } from '../../../mol-task';
import { TokenColumnProvider as TokenColumn } from '../common/text/column/token';
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
        readonly type_symbol: Column<string>,
        readonly formal_charge: Column<number>
    },
    readonly bonds: {
        readonly count: number
        readonly atomIdxA: Column<number>,
        readonly atomIdxB: Column<number>,
        readonly order: Column<number>
    }
    readonly formalCharges: {
        readonly atomIdx: Column<number>;
        readonly charge: Column<number>;
    }
}

/*
    The atom lines in a .mol file have the following structure:

    xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
    ---------------------------------------------------------------------

    Below is a breakdown of each component and its start/end indices:

    xxxxx.xxxx  (X COORDINATE, 1-10)
    yyyyy.yyyy  (Y COORDINATE, 10-20)
    zzzzz.zzzz  (Z COORDINATE, 20-30)
    _           (30 IS EMPTY)
    aaa         (ATOM SYMBOL, 31-34)
    dd          (MASS DIFF, 34-36)
    ccc         (FORMAL CHARGE, 36-39)
    sss         (ATOM STEREO PARITY, 39-42)
    hhh         (HYDROGEN COUNT+1, 42-45)
    bbb         (STEREO CARE BOX, 45-48)
    vvv         (VALENCE, 48-51)
    HHH         (H0 DESIGNATOR, 51-54)
    rrr         (UNUSED, 54-57)
    iii         (UNUSED, 57-60)
    mmm         (ATOM-ATOM MAPPING NUMBER, 60-63)
    nnn         (INVERSION/RETENTION FLAG, 63-66)
    eee         (EXACT CHANGE FLAG, 66-69)
*/
export function handleAtoms(tokenizer: Tokenizer, count: number): MolFile['atoms'] {
    const x = TokenBuilder.create(tokenizer.data, count * 2);
    const y = TokenBuilder.create(tokenizer.data, count * 2);
    const z = TokenBuilder.create(tokenizer.data, count * 2);
    const type_symbol = TokenBuilder.create(tokenizer.data, count * 2);
    const formal_charge = TokenBuilder.create(tokenizer.data, count * 2);

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
        Tokenizer.trim(tokenizer, s + 36, s + 39);
        TokenBuilder.addUnchecked(formal_charge, tokenizer.tokenStart, tokenizer.tokenEnd);
        tokenizer.position = position;
    }

    return {
        count,
        x: TokenColumn(x)(Column.Schema.float),
        y: TokenColumn(y)(Column.Schema.float),
        z: TokenColumn(z)(Column.Schema.float),
        type_symbol: TokenColumn(type_symbol)(Column.Schema.str),
        formal_charge: TokenColumn(formal_charge)(Column.Schema.int)
    };
}

export function handleBonds(tokenizer: Tokenizer, count: number): MolFile['bonds'] {
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


export function handleFormalCharges(tokenizer: Tokenizer, count: number): MolFile['formalCharges'] {
    const atomIdx = TokenBuilder.create(tokenizer.data, count * 2);
    const charge = TokenBuilder.create(tokenizer.data, count * 2);

    let i = 0;
    while (i < 100) {

        /* An attempt to explain what happens.
            Once handleFormalCharges() is called, the atom and bond sections have
            been parsed. We are now inside the properties block of the file.

            Therefore, the "pointer" of the reader is at position 0:
            M  CHG  1   2  -1
            ^
            Read the property type (positions 3 to 5):
            M  CHG  1   2  -1
            ___^^^

            If it's a charge property (CHG) we'll add it to the list of
            formal charges.
            We read the characters at positions 12 to 14 (2__),
            cleanup the spaces/tabs (2) and assign it to atomIdx property of
            the "MolFile" object.
            Same for the next triplet at positions 15 to 17.
            (-1_) becomes (-1) and is assigned to
            charge property of the "MolFile" object.
        */

        Tokenizer.markLine(tokenizer);
        const { tokenStart: s } = tokenizer;

        Tokenizer.trim(tokenizer, s + 3, s + 6);
        const propertyType = Tokenizer.getTokenString(tokenizer);

        if (propertyType === 'CHG') {
            Tokenizer.trim(tokenizer, s + 12, s + 15);
            TokenBuilder.addUnchecked(atomIdx, tokenizer.tokenStart, tokenizer.tokenEnd);
            Tokenizer.trim(tokenizer, s + 15, s + 18);
            TokenBuilder.addUnchecked(charge, tokenizer.tokenStart, tokenizer.tokenEnd);
        }
        if (propertyType === 'END') break;
        i++;
    }

    return {
        atomIdx: TokenColumn(atomIdx)(Column.Schema.int),
        charge: TokenColumn(charge)(Column.Schema.int),
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

    const formalCharges = handleFormalCharges(tokenizer, atomCount);

    const result: MolFile = {
        title,
        program,
        comment,
        atoms,
        bonds,
        formalCharges,
    };
    return Result.success(result);
}

export function parseMol(data: string) {
    return Task.create<Result<MolFile>>('Parse Mol', async () => {
        return parseInternal(data);
    });
}