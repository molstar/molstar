/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

/**
 * @param key - The value found at the atom block.
 * @returns The actual formal charge based on the mapping.
 */
export function formalChargeMapper(key: number) {
    switch (key) {
        case 7: return -3;
        case 6: return -2;
        case 5: return -1;
        case 0: return 0;
        case 3: return 1;
        case 2: return 2;
        case 1: return 3;
        case 4: return 0;
        default:
            console.error(`Value ${key} is outside the 0-7 range, defaulting to 0.`);
            return 0;
    }
}

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

interface FormalChargesRawData {
    atomIdx: Array<number>;
    charge: Array<number>;
}
export function handleFormalCharges(tokenizer: Tokenizer, lineStart: number, formalCharges: FormalChargesRawData) {

    Tokenizer.trim(tokenizer, lineStart + 6, lineStart + 9);
    const numOfCharges = parseInt(Tokenizer.getTokenString(tokenizer));
    for (let i = 0; i < numOfCharges; ++i) {
        /*
        M  CHG  3   1  -1   2   0   2  -1
                |   |   |   |   |
                |   |   |   |   |__charge2 (etc.)
                |   |   |   |
                |   |   |   |__atomIdx2
                |   |   |
                |   |   |__charge1
                |   |
                |   |__atomIdx1 (cursor at position 12)
                |
                |___numOfCharges
        */
        const offset = 9 + (i * 8);

        Tokenizer.trim(tokenizer, lineStart + offset, lineStart + offset + 4);
        const _atomIdx = Tokenizer.getTokenString(tokenizer);
        formalCharges.atomIdx.push(+_atomIdx);
        Tokenizer.trim(tokenizer, lineStart + offset + 4, lineStart + offset + 8);
        const _charge = Tokenizer.getTokenString(tokenizer);
        formalCharges.charge.push(+_charge);
    }
    /* Once the line is read, move to the next one. */
    Tokenizer.eatLine(tokenizer);
}

/** Call an appropriate handler based on the property type.
 * (For now it only calls the formal charge handler, additional handlers can
 * be added for other properties.)
 */
export function handlePropertiesBlock(tokenizer: Tokenizer): MolFile['formalCharges'] {

    const _atomIdx: Array<number> = [];
    const _charge: Array<number> = [];
    const _formalCharges: FormalChargesRawData = { atomIdx: _atomIdx, charge: _charge };

    while (tokenizer.position < tokenizer.length) {
        const { position: s } = tokenizer;

        Tokenizer.trim(tokenizer, s + 3, s + 6);
        const propertyType = Tokenizer.getTokenString(tokenizer);

        if (propertyType === 'END') break;
        Tokenizer.eatLine(tokenizer);

        switch (propertyType) {
            case 'CHG':
                handleFormalCharges(tokenizer, s, _formalCharges);
                break;
            default:
                break;
        }
    }

    const formalCharges: MolFile['formalCharges'] = {
        atomIdx: Column.ofIntArray(_formalCharges.atomIdx),
        charge: Column.ofIntArray(_formalCharges.charge)
    };
    return formalCharges;
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

    const formalCharges = handlePropertiesBlock(tokenizer);

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