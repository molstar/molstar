/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Jason Pattle <jpattle@exscientia.co.uk>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 */

import { Column } from '../../../mol-data/db';
import { MolFile, handleAtoms, handleBonds, handlePropertiesBlock } from '../mol/parser';
import { Task } from '../../../mol-task';
import { ReaderResult as Result } from '../result';
import { Tokenizer, TokenBuilder } from '../common/text/tokenizer';
import { TokenColumnProvider as TokenColumn } from '../common/text/column/token';
import { handleAtomsV3, handleBondsV3, handleCountsV3, isV3 } from './parser-v3-util';
import { StringLike } from '../../common/string-like';


/** http://c4.cabrillo.edu/404/ctfile.pdf - page 41 & 79 */

export interface SdfFileCompound {
    readonly molFile: MolFile,
    readonly dataItems: {
        readonly dataHeader: Column<string>,
        readonly data: Column<string>
    }
}

export interface SdfFile {
    readonly compounds: SdfFileCompound[]
}


const delimiter = '$$$$';

function handleDataItems(tokenizer: Tokenizer): { dataHeader: Column<string>, data: Column<string> } {
    const dataHeader = TokenBuilder.create(tokenizer.data, 32);
    const data = TokenBuilder.create(tokenizer.data, 32);

    while (tokenizer.position < tokenizer.length) {
        const line = Tokenizer.readLine(tokenizer);
        if (line.startsWith(delimiter)) break;
        if (!line) continue;

        if (line.startsWith('> ')) {
            TokenBuilder.add(dataHeader, tokenizer.tokenStart + 2, tokenizer.tokenEnd);

            Tokenizer.markLine(tokenizer);
            const start = tokenizer.tokenStart;
            let end = tokenizer.tokenEnd;
            let added = false;
            while (tokenizer.position < tokenizer.length) {
                const line2 = Tokenizer.readLine(tokenizer);
                if (!line2 || line2.startsWith(delimiter) || line2.startsWith('> ')) {
                    TokenBuilder.add(data, start, end);
                    added = true;
                    break;
                }
                end = tokenizer.tokenEnd;
            }

            if (!added) {
                TokenBuilder.add(data, start, end);
            }
        }
    }

    return {
        dataHeader: TokenColumn(dataHeader)(Column.Schema.str),
        data: TokenColumn(data)(Column.Schema.str)
    };
}

function handleCountsV2(countsAndVersion: string): { atomCount: number, bondCount: number } {
    return {
        atomCount: +countsAndVersion.substr(0, 3),
        bondCount: +countsAndVersion.substr(3, 3)
    };
}

function handleMolFile(tokenizer: Tokenizer) {
    const title = Tokenizer.readLine(tokenizer).trim();
    const program = Tokenizer.readLine(tokenizer).trim();
    const comment = Tokenizer.readLine(tokenizer).trim();

    const countsAndVersion = Tokenizer.readLine(tokenizer);
    const molIsV3 = isV3(countsAndVersion);

    const { atomCount, bondCount } = molIsV3 ? handleCountsV3(tokenizer) : handleCountsV2(countsAndVersion);

    if (Number.isNaN(atomCount) || Number.isNaN(bondCount)) {
        // try to skip to next molecule
        while (tokenizer.position < tokenizer.length) {
            const line = Tokenizer.readLine(tokenizer);
            if (line.startsWith(delimiter)) break;
        }
        return;
    }

    /* No support for formal charge parsing in V3000 molfiles at the moment,
    so all charges default to 0.*/
    const nullFormalCharges: MolFile['formalCharges'] = {
        atomIdx: Column.ofConst(0, atomCount, Column.Schema.int),
        charge: Column.ofConst(0, atomCount, Column.Schema.int)
    };

    const atoms = molIsV3 ? handleAtomsV3(tokenizer, atomCount) : handleAtoms(tokenizer, atomCount);
    const bonds = molIsV3 ? handleBondsV3(tokenizer, bondCount) : handleBonds(tokenizer, bondCount);
    const properties = molIsV3 ? { formalCharges: nullFormalCharges } : handlePropertiesBlock(tokenizer);
    const dataItems = handleDataItems(tokenizer);

    return {
        molFile: { title, program, comment, atoms, bonds, ...properties },
        dataItems
    };
}

function parseInternal(data: StringLike): Result<SdfFile> {
    const tokenizer = Tokenizer(data);

    const compounds: SdfFile['compounds'] = [];
    while (tokenizer.position < tokenizer.length) {
        const c = handleMolFile(tokenizer);
        if (c) compounds.push(c);
    }

    return Result.success({ compounds });
}

export function parseSdf(data: StringLike) {
    return Task.create<Result<SdfFile>>('Parse Sdf', async () => {
        return parseInternal(data);
    });
}