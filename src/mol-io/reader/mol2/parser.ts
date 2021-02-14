/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Zepei Xu <xuzepei19950617@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

//               NOTES
// When want to created undefined string column, must use
// undefStr = UndefinedColumn(molecule.num_atoms, ColumnType.str)
// but not
// const undefPooledStr = UndefinedColumn(molecule.num_atoms, ColumnType.pooledStr);
// because latter actuall return a column of zeros
import { Column } from '../../../mol-data/db';
import { TokenBuilder, Tokenizer } from '../common/text/tokenizer';
import { TokenColumnProvider as TokenColumn } from '../common/text/column/token';
import * as Schema from './schema';
import { ReaderResult as Result } from '../result';
import { Task, RuntimeContext, chunkedSubtask } from '../../../mol-task';

const { skipWhitespace, eatValue, markLine, getTokenString, readLine } = Tokenizer;

interface State {
    tokenizer: Tokenizer,
    molecule: Schema.Mol2Molecule,
    runtimeCtx: RuntimeContext
}

function createEmptyMolecule(): Schema.Mol2Molecule {
    return {
        mol_name: '',
        num_atoms: 0,
        num_bonds: 0,
        num_subst: 0,
        num_feat: 0,
        num_sets: 0,
        mol_type: '',
        charge_type: '',
        status_bits: '',
        mol_comment: ''
    };
}

function State(tokenizer: Tokenizer, runtimeCtx: RuntimeContext): State {
    return {
        tokenizer,
        molecule: createEmptyMolecule(),
        runtimeCtx
    };
}

const reWhitespace = /\s+/g;

function handleMolecule(state: State) {
    const { tokenizer, molecule } = state;

    while (getTokenString(tokenizer) !== '@<TRIPOS>MOLECULE' && tokenizer.position < tokenizer.data.length) {
        markLine(tokenizer);
    }

    markLine(tokenizer);
    molecule.mol_name = getTokenString(tokenizer);

    markLine(tokenizer);
    const values = getTokenString(tokenizer).trim().split(reWhitespace);
    molecule.num_atoms = parseInt(values[0]);
    molecule.num_bonds = parseInt(values[1]);
    molecule.num_subst = parseInt(values[2]);
    molecule.num_feat = parseInt(values[3]);
    molecule.num_sets = parseInt(values[4]);

    markLine(tokenizer);
    const mol_type = getTokenString(tokenizer);
    if (mol_type.startsWith('@<TRIPOS>')) return;
    molecule.mol_type = mol_type;

    markLine(tokenizer);
    const charge_type = getTokenString(tokenizer);
    if (charge_type.startsWith('@<TRIPOS>')) return;
    molecule.charge_type = charge_type;

    markLine(tokenizer);
    const status_bits = getTokenString(tokenizer);
    if (status_bits.startsWith('@<TRIPOS>')) return;
    molecule.status_bits = status_bits;

    markLine(tokenizer);
    const mol_comment = getTokenString(tokenizer);
    if (mol_comment.startsWith('@<TRIPOS>')) return;
    molecule.mol_comment = mol_comment;
}

async function handleAtoms(state: State): Promise<Schema.Mol2Atoms> {
    const { tokenizer, molecule } = state;

    // skip empty lines and '@<TRIPOS>ATOM'
    while (getTokenString(tokenizer) !== '@<TRIPOS>ATOM' && tokenizer.position < tokenizer.data.length) {
        markLine(tokenizer);
    }

    const initialTokenizerPosition = tokenizer.position;
    const initialTokenizerLineNumber = tokenizer.lineNumber;
    const firstLine = readLine(tokenizer);
    const firstLineArray = firstLine.trim().split(/\s+/g);
    const columnCount = firstLineArray.length;

    // columns
    const atom_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const atom_nameTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const xTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const yTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const zTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const atom_typeTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const subst_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const subst_nameTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const chargeTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const status_bitTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);

    const undefFloat = Column.Undefined(molecule.num_atoms, Column.Schema.float);
    const undefInt = Column.Undefined(molecule.num_atoms, Column.Schema.int);
    const undefStr = Column.Undefined(molecule.num_atoms, Column.Schema.str);

    tokenizer.position = initialTokenizerPosition;
    tokenizer.lineNumber = initialTokenizerLineNumber;

    const { length } = tokenizer;
    let linesAlreadyRead = 0;
    await chunkedSubtask(state.runtimeCtx, 100000, void 0, chunkSize => {
        const linesToRead = Math.min(molecule.num_atoms - linesAlreadyRead, chunkSize);
        for (let i = 0; i < linesToRead; i++) {
            for (let j = 0; j < columnCount; j++) {
                skipWhitespace(tokenizer);
                tokenizer.tokenStart = tokenizer.position;
                eatValue(tokenizer);
                switch (j) {
                    case 0:
                        TokenBuilder.addUnchecked(atom_idTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 1:
                        TokenBuilder.addUnchecked(atom_nameTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 2:
                        TokenBuilder.addUnchecked(xTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 3:
                        TokenBuilder.addUnchecked(yTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 4:
                        TokenBuilder.addUnchecked(zTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 5:
                        TokenBuilder.addUnchecked(atom_typeTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 6:
                        TokenBuilder.addUnchecked(subst_idTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 7:
                        TokenBuilder.addUnchecked(subst_nameTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 8:
                        TokenBuilder.addUnchecked(chargeTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 9:
                        TokenBuilder.addUnchecked(status_bitTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                }
            }
        }
        linesAlreadyRead += linesToRead;
        return linesToRead;
    }, ctx => ctx.update({ message: 'Parsing...', current: tokenizer.position, max: length }));

    const ret = {
        count: molecule.num_atoms,

        atom_id: TokenColumn(atom_idTokens)(Column.Schema.int),
        atom_name: TokenColumn(atom_nameTokens)(Column.Schema.str),
        x: TokenColumn(xTokens)(Column.Schema.float),
        y: TokenColumn(yTokens)(Column.Schema.float),
        z: TokenColumn(zTokens)(Column.Schema.float),

        atom_type: columnCount > 5 ? TokenColumn(atom_typeTokens)(Column.Schema.str) : undefStr,
        subst_id: columnCount > 6 ? TokenColumn(subst_idTokens)(Column.Schema.int) : undefInt,
        subst_name: columnCount > 7 ? TokenColumn(subst_nameTokens)(Column.Schema.str) : undefStr,
        charge: columnCount > 8 ? TokenColumn(chargeTokens)(Column.Schema.float) : undefFloat,
        status_bit: columnCount > 9 ? TokenColumn(status_bitTokens)(Column.Schema.str) : undefStr,
    };
    return ret;
}

async function handleBonds(state: State): Promise<Schema.Mol2Bonds> {
    const { tokenizer, molecule } = state;

    while (getTokenString(tokenizer) !== '@<TRIPOS>BOND' && tokenizer.position < tokenizer.data.length) {
        markLine(tokenizer);
    }

    const initialTokenizerPosition = tokenizer.position;
    const initialTokenizerLineNumber = tokenizer.lineNumber;
    const firstLine = readLine(tokenizer);
    const firstLineArray = firstLine.trim().split(/\s+/g);
    const columnCount = firstLineArray.length;

    // columns
    const bond_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const origin_bond_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const target_bond_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const bondTypeTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const status_bitTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);

    tokenizer.position = initialTokenizerPosition;
    tokenizer.lineNumber = initialTokenizerLineNumber;

    const { length } = tokenizer;
    let linesAlreadyRead = 0;
    await chunkedSubtask(state.runtimeCtx, 100000, void 0, chunkSize => {
        const linesToRead = Math.min(molecule.num_bonds - linesAlreadyRead, chunkSize);
        for (let i = 0; i < linesToRead; i++) {
            for (let j = 0; j < columnCount; j++) {
                skipWhitespace(tokenizer);
                tokenizer.tokenStart = tokenizer.position;
                eatValue(tokenizer);
                switch (j) {
                    case 0:
                        TokenBuilder.addUnchecked(bond_idTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 1:
                        TokenBuilder.addUnchecked(origin_bond_idTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 2:
                        TokenBuilder.addUnchecked(target_bond_idTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    case 3:
                        TokenBuilder.addUnchecked(bondTypeTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                    default:
                        TokenBuilder.addUnchecked(status_bitTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                        break;
                }
            }
        }
        linesAlreadyRead += linesToRead;
        return linesToRead;
    }, ctx => ctx.update({ message: 'Parsing...', current: tokenizer.position, max: length }));

    const ret = {
        count: molecule.num_bonds,

        bond_id: TokenColumn(bond_idTokens)(Column.Schema.int),
        origin_atom_id: TokenColumn(origin_bond_idTokens)(Column.Schema.int),
        target_atom_id: TokenColumn(target_bond_idTokens)(Column.Schema.int),
        bond_type: TokenColumn(bondTypeTokens)(Column.Schema.str),

        status_bits: columnCount > 4
            ? TokenColumn(status_bitTokens)(Column.Schema.str)
            : Column.Undefined(molecule.num_bonds, Column.Schema.str),
    };

    return ret;
}

async function parseInternal(ctx: RuntimeContext, data: string, name: string): Promise<Result<Schema.Mol2File>> {
    const tokenizer = Tokenizer(data);

    ctx.update({ message: 'Parsing...', current: 0, max: data.length });
    const structures: Schema.Mol2Structure[] = [];
    while (tokenizer.position < data.length) {
        const state = State(tokenizer, ctx);
        handleMolecule(state);
        const atoms = await handleAtoms(state);
        const bonds = await handleBonds(state);
        structures.push({ molecule: state.molecule, atoms, bonds });
        skipWhitespace(tokenizer);
        while (getTokenString(tokenizer) !== '@<TRIPOS>MOLECULE' && tokenizer.position < tokenizer.data.length) {
            markLine(tokenizer);
        }
    }

    const result: Schema.Mol2File = { name, structures };
    return Result.success(result);
}

export function parseMol2(data: string, name: string) {
    return Task.create<Result<Schema.Mol2File>>('Parse MOL2', async ctx => {
        return await parseInternal(ctx, data, name);
    });
}
