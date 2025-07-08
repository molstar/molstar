/**
 * Copyright (c) 2017-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Zepei Xu <xuzepei19950617@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Eric E <etongfu@@outlook.com>
 */

import { Column } from '../../../mol-data/db';
import { TokenBuilder, Tokenizer, Tokens } from '../common/text/tokenizer';
import { TokenColumnProvider as TokenColumn } from '../common/text/column/token';
import { ReaderResult as Result } from '../result';
import { Task, RuntimeContext, chunkedSubtask } from '../../../mol-task';
import { StringLike } from '../../common/string-like';
import { Mol2Atoms, Mol2Bonds, Mol2File, Mol2Molecule, Mol2Structure, Mol2Substructure } from './schema';

const { skipWhitespace, eatValue, markLine, getTokenString, skipStrictWhitespace } = Tokenizer;

interface State {
    tokenizer: Tokenizer,
    molecule: Mol2Molecule,
    runtimeCtx: RuntimeContext
}

function createEmptyMolecule(): Mol2Molecule {
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

/**
 * Just read the columns and get the max count of columns for the **atoms** and **bonds**
 * @param linesToRead The total lines
 * @param tokenIndexColums
 * @param tokenizer
 * @description
 * !!!This function has side effects!!!
 * 1. Called inside of the `chunkedSubtask`
 * 2. It will change the parameters of the `tokenIndexColums` , `tokenizer` and `linesToRead`
 * @returns The max count of columns
 */
function _readColumnsAndGetMaxCount(linesToRead: number, tokenIndexColums: Tokens[], tokenizer: Tokenizer): number {
    let maxColumnCount = 0;
    for (let i = 0; i < linesToRead; i++) {
        let tokenIndex = 0;
        skipWhitespace(tokenizer);
        while (true) {
            skipStrictWhitespace(tokenizer);
            tokenizer.tokenStart = tokenizer.position;
            eatValue(tokenizer);
            if (tokenizer.tokenStart === tokenizer.tokenEnd) break;
            const col = tokenIndexColums[tokenIndex++];
            if (!col) continue;
            TokenBuilder.addUnchecked(col, tokenizer.tokenStart, tokenizer.tokenEnd);
        }
        if (tokenIndex > maxColumnCount) maxColumnCount = tokenIndex;
        for (let cI = tokenIndex; cI < tokenIndexColums.length; cI++) {
            TokenBuilder.addUnchecked(tokenIndexColums[cI], 0, 0);
        }
    }
    return maxColumnCount + 1;
}

async function handleAtoms(state: State): Promise<Mol2Atoms> {
    const { tokenizer, molecule } = state;

    // skip empty lines and '@<TRIPOS>ATOM'
    while (getTokenString(tokenizer) !== '@<TRIPOS>ATOM' && tokenizer.position < tokenizer.data.length) {
        markLine(tokenizer);
    }

    const initialTokenizerPosition = tokenizer.position;
    const initialTokenizerLineNumber = tokenizer.lineNumber;

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
    const status_bitsTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);

    const undefFloat = Column.Undefined(molecule.num_atoms, Column.Schema.float);
    const undefInt = Column.Undefined(molecule.num_atoms, Column.Schema.int);
    const undefStr = Column.Undefined(molecule.num_atoms, Column.Schema.str);

    tokenizer.position = initialTokenizerPosition;
    tokenizer.lineNumber = initialTokenizerLineNumber;

    let maxColumnCount = 0;
    const tokenIndexToColumn = [
        atom_idTokens,
        atom_nameTokens,
        xTokens,
        yTokens,
        zTokens,
        atom_typeTokens,
        subst_idTokens,
        subst_nameTokens,
        chargeTokens,
        status_bitsTokens
    ];

    const { length } = tokenizer;
    let linesAlreadyRead = 0;

    await chunkedSubtask(state.runtimeCtx, 100000, void 0, chunkSize => {
        const linesToRead = Math.min(molecule.num_atoms - linesAlreadyRead, chunkSize);
        maxColumnCount = Math.max(maxColumnCount, _readColumnsAndGetMaxCount(linesToRead, tokenIndexToColumn, tokenizer));
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

        atom_type: maxColumnCount > 5 ? TokenColumn(atom_typeTokens)(Column.Schema.str) : undefStr,
        subst_id: maxColumnCount > 6 ? TokenColumn(subst_idTokens)(Column.Schema.int) : undefInt,
        subst_name: maxColumnCount > 7 ? TokenColumn(subst_nameTokens)(Column.Schema.str) : undefStr,
        charge: maxColumnCount > 8 ? TokenColumn(chargeTokens)(Column.Schema.float) : undefFloat,
        status_bits: maxColumnCount > 9 ? TokenColumn(status_bitsTokens)(Column.Schema.str) : undefStr,
    };
    return ret;
}

async function handleBonds(state: State): Promise<Mol2Bonds> {
    const { tokenizer, molecule } = state;

    while (getTokenString(tokenizer) !== '@<TRIPOS>BOND' && tokenizer.position < tokenizer.data.length) {
        markLine(tokenizer);
    }

    const initialTokenizerPosition = tokenizer.position;
    const initialTokenizerLineNumber = tokenizer.lineNumber;

    // columns
    const bond_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const origin_bond_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const target_bond_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const bondTypeTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const status_bitsTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);

    const undefStr = Column.Undefined(molecule.num_bonds, Column.Schema.str);

    tokenizer.position = initialTokenizerPosition;
    tokenizer.lineNumber = initialTokenizerLineNumber;

    let maxColumnCount = 0;
    const tokenIndexToColumn = [
        bond_idTokens,
        origin_bond_idTokens,
        target_bond_idTokens,
        bondTypeTokens,
        status_bitsTokens
    ];

    const { length } = tokenizer;
    let linesAlreadyRead = 0;

    await chunkedSubtask(state.runtimeCtx, 100000, void 0, chunkSize => {
        const linesToRead = Math.min(molecule.num_bonds - linesAlreadyRead, chunkSize);
        maxColumnCount = Math.max(maxColumnCount, _readColumnsAndGetMaxCount(linesToRead, tokenIndexToColumn, tokenizer));
        linesAlreadyRead += linesToRead;
        return linesToRead;
    }, ctx => ctx.update({ message: 'Parsing...', current: tokenizer.position, max: length }));

    const ret = {
        count: molecule.num_bonds,

        bond_id: TokenColumn(bond_idTokens)(Column.Schema.int),
        origin_atom_id: TokenColumn(origin_bond_idTokens)(Column.Schema.int),
        target_atom_id: TokenColumn(target_bond_idTokens)(Column.Schema.int),
        bond_type: TokenColumn(bondTypeTokens)(Column.Schema.str),

        status_bits: maxColumnCount > 4 ? TokenColumn(status_bitsTokens)(Column.Schema.str) : undefStr,
    };

    return ret;
}

async function handleSubstructures(state: State): Promise<Mol2Substructure | undefined> {
    const { tokenizer, molecule } = state;
    if (molecule.num_subst === 0) return;

    while (getTokenString(tokenizer) !== '@<TRIPOS>SUBSTRUCTURE' && tokenizer.position < tokenizer.data.length) {
        markLine(tokenizer);
    }

    const initialTokenizerPosition = tokenizer.position;
    const initialTokenizerLineNumber = tokenizer.lineNumber;

    // columns
    const subst_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_subst * 2);
    const subst_nameTokens = TokenBuilder.create(tokenizer.data, molecule.num_subst * 2);
    const root_atomTokens = TokenBuilder.create(tokenizer.data, molecule.num_subst * 2);
    const subst_typeTokens = TokenBuilder.create(tokenizer.data, molecule.num_subst * 2);
    const dict_typeTokens = TokenBuilder.create(tokenizer.data, molecule.num_subst * 2);
    const chainTokens = TokenBuilder.create(tokenizer.data, molecule.num_subst * 2);
    const sub_typeTokens = TokenBuilder.create(tokenizer.data, molecule.num_subst * 2);
    const inter_bondsTokens = TokenBuilder.create(tokenizer.data, molecule.num_subst * 2);
    const status_bitsTokens = TokenBuilder.create(tokenizer.data, molecule.num_subst * 2);

    const undefInt = Column.Undefined(molecule.num_subst, Column.Schema.int);
    const undefStr = Column.Undefined(molecule.num_subst, Column.Schema.str);

    tokenizer.position = initialTokenizerPosition;
    tokenizer.lineNumber = initialTokenizerLineNumber;

    let maxColumnCount = 0;
    const tokenIndexToColumn = [
        subst_idTokens,
        subst_nameTokens,
        root_atomTokens,
        subst_typeTokens,
        dict_typeTokens,
        chainTokens,
        sub_typeTokens,
        inter_bondsTokens,
        status_bitsTokens
    ];

    const { length } = tokenizer;
    let linesAlreadyRead = 0;

    await chunkedSubtask(state.runtimeCtx, 100000, void 0, chunkSize => {
        const linesToRead = Math.min(molecule.num_subst - linesAlreadyRead, chunkSize);
        maxColumnCount = Math.max(maxColumnCount, _readColumnsAndGetMaxCount(linesToRead, tokenIndexToColumn, tokenizer));
        linesAlreadyRead += linesToRead;
        return linesToRead;
    }, ctx => ctx.update({ message: 'Parsing...', current: tokenizer.position, max: length }));

    const ret = {
        count: molecule.num_subst,

        subst_id: TokenColumn(subst_idTokens)(Column.Schema.int),
        subst_name: TokenColumn(subst_nameTokens)(Column.Schema.str),
        root_atom: TokenColumn(root_atomTokens)(Column.Schema.int),

        subst_type: maxColumnCount > 3 ? TokenColumn(subst_typeTokens)(Column.Schema.str) : undefStr,
        dict_type: maxColumnCount > 4 ? TokenColumn(dict_typeTokens)(Column.Schema.str) : undefStr,
        chain: maxColumnCount > 5 ? TokenColumn(chainTokens)(Column.Schema.str) : undefStr,
        sub_type: maxColumnCount > 6 ? TokenColumn(sub_typeTokens)(Column.Schema.str) : undefStr,
        inter_bonds: maxColumnCount > 7 ? TokenColumn(inter_bondsTokens)(Column.Schema.int) : undefInt,
        status_bits: maxColumnCount > 8 ? TokenColumn(status_bitsTokens)(Column.Schema.str) : undefStr,
    };

    return ret;
}

function handleCrysin(state: State) {
    const { tokenizer } = state;

    while (tokenizer.position < tokenizer.data.length) {
        const l = getTokenString(tokenizer);
        if (l === '@<TRIPOS>MOLECULE') {
            return;
        } else if (l === '@<TRIPOS>CRYSIN') {
            break;
        } else {
            markLine(tokenizer);
        }
    }

    if (tokenizer.position >= tokenizer.data.length) return;

    markLine(tokenizer);
    const values = getTokenString(tokenizer).trim().split(reWhitespace);
    return {
        a: parseFloat(values[0]),
        b: parseFloat(values[1]),
        c: parseFloat(values[2]),
        alpha: parseFloat(values[3]),
        beta: parseFloat(values[4]),
        gamma: parseFloat(values[5]),
        spaceGroup: parseInt(values[6], 10),
        setting: parseInt(values[7], 10),
    };
}

async function parseInternal(ctx: RuntimeContext, data: StringLike, name: string): Promise<Result<Mol2File>> {
    const tokenizer = Tokenizer(data);

    ctx.update({ message: 'Parsing...', current: 0, max: data.length });
    const structures: Mol2Structure[] = [];
    while (tokenizer.position < data.length) {
        const state = State(tokenizer, ctx);
        handleMolecule(state);
        const atoms = await handleAtoms(state);
        const bonds = await handleBonds(state);
        const substructures = await handleSubstructures(state);
        const crysin = handleCrysin(state);
        structures.push({ molecule: state.molecule, atoms, bonds, substructures, crysin });
        skipWhitespace(tokenizer);
        while (getTokenString(tokenizer) !== '@<TRIPOS>MOLECULE' && tokenizer.position < tokenizer.data.length) {
            markLine(tokenizer);
        }
    }

    const result: Mol2File = { name, structures };
    return Result.success(result);
}

export function parseMol2(data: StringLike, name: string) {
    return Task.create<Result<Mol2File>>('Parse MOL2', async ctx => {
        return await parseInternal(ctx, data, name);
    });
}
