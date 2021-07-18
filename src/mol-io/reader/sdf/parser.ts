/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '../../../mol-data/db';
import { MolFile, handleAtoms, handleBonds } from '../mol/parser';
import { Task } from '../../../mol-task';
import { ReaderResult as Result } from '../result';
import { Tokenizer, TokenBuilder } from '../common/text/tokenizer';
import { TokenColumnProvider as TokenColumn } from '../common/text/column/token';

/** http://c4.cabrillo.edu/404/ctfile.pdf - page 41 */

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

function handleMolFile(tokenizer: Tokenizer) {
    const title = Tokenizer.readLine(tokenizer).trim();
    const program = Tokenizer.readLine(tokenizer).trim();
    const comment = Tokenizer.readLine(tokenizer).trim();

    const counts = Tokenizer.readLine(tokenizer);

    const atomCount = +counts.substr(0, 3), bondCount = +counts.substr(3, 3);

    if (Number.isNaN(atomCount) || Number.isNaN(bondCount)) {
        // try to skip to next molecule
        while (tokenizer.position < tokenizer.length) {
            const line = Tokenizer.readLine(tokenizer);
            if (line.startsWith(delimiter)) break;
        }
        return;
    }

    const atoms = handleAtoms(tokenizer, atomCount);
    const bonds = handleBonds(tokenizer, bondCount);
    const dataItems = handleDataItems(tokenizer);

    return {
        molFile: { title, program, comment, atoms, bonds },
        dataItems
    };
}

function parseInternal(data: string): Result<SdfFile> {
    const tokenizer = Tokenizer(data);

    const compounds: SdfFile['compounds'] = [];
    while (tokenizer.position < tokenizer.length) {
        const c = handleMolFile(tokenizer);
        if (c) compounds.push(c);
    }

    return Result.success({ compounds });
}

export function parseSdf(data: string) {
    return Task.create<Result<SdfFile>>('Parse Sdf', async () => {
        return parseInternal(data);
    });
}