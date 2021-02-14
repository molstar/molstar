/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Column } from '../../../mol-data/db';
import { MolFile, handleAtoms, handleBonds } from '../mol/parser';
import { Task } from '../../../mol-task';
import { ReaderResult as Result } from '../result';
import { Tokenizer, TokenBuilder } from '../common/text/tokenizer';
import { TokenColumnProvider as TokenColumn } from '../common/text/column/token';

/** http://c4.cabrillo.edu/404/ctfile.pdf - page 41 */
export interface SdfFile {
    readonly compounds: {
        readonly molFile: MolFile,
        readonly dataItems: {
            readonly dataHeader: Column<string>,
            readonly data: Column<string>
        }
    }[]
}

function handleDataItems(tokenizer: Tokenizer): { dataHeader: Column<string>, data: Column<string> } {
    const dataHeader = TokenBuilder.create(tokenizer.data, 32);
    const data = TokenBuilder.create(tokenizer.data, 32);

    let sawHeaderToken = false;
    while (tokenizer.position < tokenizer.length) {
        const line = Tokenizer.readLine(tokenizer);
        if (!!line) {
            if (line.startsWith('> <')) {
                TokenBuilder.add(dataHeader, tokenizer.tokenStart + 3, tokenizer.tokenEnd - 1);
                sawHeaderToken = true;
            } else if (sawHeaderToken) {
                TokenBuilder.add(data, tokenizer.tokenStart, tokenizer.tokenEnd);
                sawHeaderToken = false;
                // TODO can there be multiline values?
            }
        } else {
            sawHeaderToken = false;
        }
    }

    return {
        dataHeader: TokenColumn(dataHeader)(Column.Schema.str),
        data: TokenColumn(data)(Column.Schema.str)
    };
}

function handleMolFile(data: string) {
    const tokenizer = Tokenizer(data);

    const title = Tokenizer.readLine(tokenizer).trim();
    const program = Tokenizer.readLine(tokenizer).trim();
    const comment = Tokenizer.readLine(tokenizer).trim();

    const counts = Tokenizer.readLine(tokenizer);

    const atomCount = +counts.substr(0, 3), bondCount = +counts.substr(3, 3);

    const atoms = handleAtoms(tokenizer, atomCount);
    const bonds = handleBonds(tokenizer, bondCount);
    const dataItems = handleDataItems(tokenizer);

    return {
        molFile: { title, program, comment, atoms, bonds },
        dataItems
    };
}

const delimiter = '$$$$';
function parseInternal(data: string): Result<SdfFile> {
    const result: SdfFile = { compounds: data.split(delimiter).map(d => handleMolFile(d)) };
    return Result.success(result);
}

export function parseSdf(data: string) {
    return Task.create<Result<SdfFile>>('Parse Sdf', async () => {
        return parseInternal(data);
    });
}