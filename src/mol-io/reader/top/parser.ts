/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from '../../../mol-task';
import { Tokenizer, TokenBuilder } from '../common/text/tokenizer';
import { ReaderResult as Result } from '../result';
import { TokenColumnProvider as TokenColumn } from '../common/text/column/token';
import { Column, Table } from '../../../mol-data/db';
import { Mutable } from '../../../mol-util/type-helpers';
import { StringLike } from '../../common/string-like';


// https://manual.gromacs.org/2021-current/reference-manual/file-formats.html#top

const AtomsSchema = {
    nr: Column.Schema.Int(),
    type: Column.Schema.Str(),
    resnr: Column.Schema.Int(),
    residu: Column.Schema.Str(),
    atom: Column.Schema.Str(),
    cgnr: Column.Schema.Int(),
    charge: Column.Schema.Float(),
    mass: Column.Schema.Float(),
};

const BondsSchema = {
    ai: Column.Schema.Int(),
    aj: Column.Schema.Int(),
};

const MoleculesSchema = {
    compound: Column.Schema.Str(),
    molCount: Column.Schema.Int(),
};

type Compound = {
    atoms: Table<typeof AtomsSchema>
    bonds?: Table<typeof BondsSchema>
}

export interface TopFile {
    readonly system: string
    readonly molecules: Table<typeof MoleculesSchema>
    readonly compounds: Record<string, Compound>
}

const { readLine, markLine, skipWhitespace, markStart, eatValue, eatLine } = Tokenizer;

function State(tokenizer: Tokenizer, runtimeCtx: RuntimeContext) {
    return {
        tokenizer,
        runtimeCtx,
    };
}
type State = ReturnType<typeof State>

const reField = /\[ (.+) \]/;
const reWhitespace = /\s+/;

function handleMoleculetype(state: State) {
    const { tokenizer } = state;

    let molName: string | undefined = undefined;

    while (tokenizer.tokenEnd < tokenizer.length) {
        skipWhitespace(tokenizer);
        const c = tokenizer.data.charAt(tokenizer.position);
        if (c === '[') break;
        if (c === ';' || c === '*') {
            markLine(tokenizer);
            continue;
        }

        if (molName !== undefined) throw new Error('more than one molName');

        const line = readLine(tokenizer);
        molName = line.split(reWhitespace)[0];
    }

    if (molName === undefined) throw new Error('missing molName');

    return molName;
}

function handleAtoms(state: State) {
    const { tokenizer } = state;

    const nr = TokenBuilder.create(tokenizer.data, 64);
    const type = TokenBuilder.create(tokenizer.data, 64);
    const resnr = TokenBuilder.create(tokenizer.data, 64);
    const residu = TokenBuilder.create(tokenizer.data, 64);
    const atom = TokenBuilder.create(tokenizer.data, 64);
    const cgnr = TokenBuilder.create(tokenizer.data, 64);
    const charge = TokenBuilder.create(tokenizer.data, 64);
    const mass = TokenBuilder.create(tokenizer.data, 64);

    while (tokenizer.tokenEnd < tokenizer.length) {
        skipWhitespace(tokenizer);
        const c = tokenizer.data.charAt(tokenizer.position);
        if (c === '[') break;
        if (c === ';' || c === '*') {
            markLine(tokenizer);
            continue;
        }

        for (let j = 0; j < 8; ++j) {
            skipWhitespace(tokenizer);
            markStart(tokenizer);
            eatValue(tokenizer);

            switch (j) {
                case 0: TokenBuilder.add(nr, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                case 1: TokenBuilder.add(type, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                case 2: TokenBuilder.add(resnr, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                case 3: TokenBuilder.add(residu, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                case 4: TokenBuilder.add(atom, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                case 5: TokenBuilder.add(cgnr, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                case 6: TokenBuilder.add(charge, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                case 7: TokenBuilder.add(mass, tokenizer.tokenStart, tokenizer.tokenEnd); break;
            }
        }
        // ignore any extra columns
        markLine(tokenizer);
    }

    return Table.ofColumns(AtomsSchema, {
        nr: TokenColumn(nr)(Column.Schema.int),
        type: TokenColumn(type)(Column.Schema.str),
        resnr: TokenColumn(resnr)(Column.Schema.int),
        residu: TokenColumn(residu)(Column.Schema.str),
        atom: TokenColumn(atom)(Column.Schema.str),
        cgnr: TokenColumn(cgnr)(Column.Schema.int),
        charge: TokenColumn(charge)(Column.Schema.float),
        mass: TokenColumn(mass)(Column.Schema.float),
    });
}

function handleBonds(state: State) {
    const { tokenizer } = state;

    const ai = TokenBuilder.create(tokenizer.data, 64);
    const aj = TokenBuilder.create(tokenizer.data, 64);

    while (tokenizer.tokenEnd < tokenizer.length) {
        skipWhitespace(tokenizer);
        const c = tokenizer.data.charAt(tokenizer.position);
        if (c === '[') break;
        if (c === ';' || c === '*') {
            markLine(tokenizer);
            continue;
        }

        for (let j = 0; j < 2; ++j) {
            skipWhitespace(tokenizer);
            markStart(tokenizer);
            eatValue(tokenizer);

            switch (j) {
                case 0: TokenBuilder.add(ai, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                case 1: TokenBuilder.add(aj, tokenizer.tokenStart, tokenizer.tokenEnd); break;
            }
        }
        // ignore any extra columns
        markLine(tokenizer);
    }

    return Table.ofColumns(BondsSchema, {
        ai: TokenColumn(ai)(Column.Schema.int),
        aj: TokenColumn(aj)(Column.Schema.int),
    });
}

function handleSystem(state: State) {
    const { tokenizer } = state;

    let system: string | undefined = undefined;

    while (tokenizer.tokenEnd < tokenizer.length) {
        skipWhitespace(tokenizer);
        const c = tokenizer.data.charAt(tokenizer.position);
        if (c === '[') break;
        if (c === ';' || c === '*') {
            markLine(tokenizer);
            continue;
        }

        if (system !== undefined) throw new Error('more than one system');
        system = readLine(tokenizer).trim();
    }

    if (system === undefined) throw new Error('missing system');

    return system;
}

function handleMolecules(state: State) {
    const { tokenizer } = state;

    const compound = TokenBuilder.create(tokenizer.data, 64);
    const molCount = TokenBuilder.create(tokenizer.data, 64);

    while (tokenizer.tokenEnd < tokenizer.length) {
        skipWhitespace(tokenizer);
        if (tokenizer.position >= tokenizer.length) break;

        const c = tokenizer.data.charAt(tokenizer.position);
        if (c === '[') break;
        if (c === ';' || c === '*') {
            markLine(tokenizer);
            continue;
        }

        for (let j = 0; j < 2; ++j) {
            skipWhitespace(tokenizer);
            markStart(tokenizer);
            eatValue(tokenizer);

            switch (j) {
                case 0: TokenBuilder.add(compound, tokenizer.tokenStart, tokenizer.tokenEnd); break;
                case 1: TokenBuilder.add(molCount, tokenizer.tokenStart, tokenizer.tokenEnd); break;
            }
        }
        // ignore any extra columns
        eatLine(tokenizer);
        markStart(tokenizer);
    }

    return Table.ofColumns(MoleculesSchema, {
        compound: TokenColumn(compound)(Column.Schema.str),
        molCount: TokenColumn(molCount)(Column.Schema.int),
    });
}

async function parseInternal(data: StringLike, ctx: RuntimeContext): Promise<Result<TopFile>> {
    const t = Tokenizer(data);
    const state = State(t, ctx);

    const result: Mutable<TopFile> = Object.create(null);
    let prevPosition = 0;

    result.compounds = {};
    let currentCompound: Partial<Compound> = {};
    let currentMolName = '';

    function addMol() {
        if (currentMolName && currentCompound.atoms) {
            result.compounds[currentMolName] = currentCompound as Compound;
            currentCompound = {};
            currentMolName = '';
        }
    }

    while (t.tokenEnd < t.length) {
        if (t.position - prevPosition > 100000 && ctx.shouldUpdate) {
            prevPosition = t.position;
            await ctx.update({ current: t.position, max: t.length });
        }

        const line = readLine(state.tokenizer).trim();

        if (!line || line[0] === '*' || line[0] === ';') {
            continue;
        }

        if (line.startsWith('#include')) {
            throw new Error('#include statements not allowed');
        }

        if (line.startsWith('[')) {
            const fieldMatch = line.match(reField);
            if (fieldMatch === null) throw new Error('expected field name');

            const fieldName = fieldMatch[1];
            if (fieldName === 'moleculetype') {
                addMol();
                currentMolName = handleMoleculetype(state);
            } else if (fieldName === 'atoms') {
                currentCompound.atoms = handleAtoms(state);
            } else if (fieldName === 'bonds') {
                currentCompound.bonds = handleBonds(state);
            } else if (fieldName === 'system') {
                result.system = handleSystem(state);
            } else if (fieldName === 'molecules') {
                addMol(); // add the last compound
                result.molecules = handleMolecules(state);
            } else {
                while (t.tokenEnd < t.length) {
                    if (t.data.charAt(t.position) === '[') break;
                    markLine(t);
                }
            }
        }
    }

    return Result.success(result);
}

export function parseTop(data: StringLike) {
    return Task.create<Result<TopFile>>('Parse TOP', async ctx => {
        return await parseInternal(data, ctx);
    });
}