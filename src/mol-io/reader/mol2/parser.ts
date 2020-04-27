/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import TokenColumn from '../common/text/column/token';
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
    molecule.mol_type = getTokenString(tokenizer);

    markLine(tokenizer);
    molecule.charge_type = getTokenString(tokenizer);

    markLine(tokenizer);
    if (getTokenString(tokenizer) === '') return;
    molecule.status_bits = getTokenString(tokenizer);

    markLine(tokenizer);
    if (getTokenString(tokenizer) === '') return;
    molecule.mol_comment = getTokenString(tokenizer);
}

function isStatus_bit(aString: string): boolean {
    if (aString.includes('DSPMOD') || aString.includes('TYPECOL') || aString.includes('CAP')
        || aString.includes('BACKBONE') || aString.includes('DICT') || aString.includes('ESSENTIAL')
        || aString.includes('WATER') || aString.includes('DIRECT')) {
        return true;
    }
    return false;
}

async function handleAtoms(state: State): Promise<Schema.Mol2Atoms> {
    const { tokenizer, molecule } = state;
    let hasSubst_id = false;
    let hasSubst_name = false;
    let hasCharge = false;
    let hasStatus_bit = false;

    // skip empty lines and '@<TRIPOS>ATOM'
    while (getTokenString(tokenizer) !== '@<TRIPOS>ATOM' && tokenizer.position < tokenizer.data.length) {
        markLine(tokenizer);
    }

    const initialTokenizerPosition = tokenizer.position;
    const initialTokenizerLineNumber = tokenizer.lineNumber;
    const firstLine = readLine(tokenizer);
    const firstLineArray = firstLine.trim().split(/\s+/g);
    const firstLineLength = firstLineArray.length;

    // optional columns are in order "integer string float string".
    // Use this to find out which column is missing or empty
    for (let i = 6; i < firstLineLength; i++) {
        if (!isNaN(Number(firstLineArray[i]))) {
            if (firstLineArray[i].indexOf('.') === -1) {
                hasSubst_id = true;
            } else {
                hasCharge = true;
            }
        } else if (isNaN(Number(firstLineArray[i]))) {
            if (!isStatus_bit(firstLineArray[i])) {
                hasSubst_name = true;
            } else {
                hasStatus_bit = true;
            }
        }
    }

    // required columns
    const atom_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const atom_nameTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const xTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const yTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const zTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const atom_typeTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);

    const atom_idTokenColumn = TokenColumn(atom_idTokens);
    const atom_nameTokenColumn = TokenColumn(atom_nameTokens);
    const xTokenColumn = TokenColumn(xTokens);
    const yTokenColumn = TokenColumn(yTokens);
    const zTokenColumn = TokenColumn(zTokens);
    const atom_typeColumn = TokenColumn(atom_typeTokens);

    // optional columns
    const subst_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const subst_nameTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const chargeTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);
    const status_bitTokens = TokenBuilder.create(tokenizer.data, molecule.num_atoms * 2);

    const subst_idTokenColumn = TokenColumn(subst_idTokens);
    const subst_nameTokenColumn = TokenColumn(subst_nameTokens);
    const chargeTokenColumn = TokenColumn(chargeTokens);
    const status_bitTokenColumn = TokenColumn(status_bitTokens);

    const undefFloat = Column.Undefined(molecule.num_atoms, Column.Schema.float);
    const undefInt = Column.Undefined(molecule.num_atoms, Column.Schema.int);
    const undefStr = Column.Undefined(molecule.num_atoms, Column.Schema.str);

    let numOfColumn = 6;
    if (hasSubst_id) { numOfColumn++; }
    if (hasSubst_name) { numOfColumn++; }
    if (hasCharge) { numOfColumn++; }
    if (hasStatus_bit) { numOfColumn++; }

    tokenizer.position = initialTokenizerPosition;
    tokenizer.lineNumber = initialTokenizerLineNumber;

    const { length } = tokenizer;
    let linesAlreadyRead = 0;
    await chunkedSubtask(state.runtimeCtx, 100000, void 0, chunkSize => {
        const linesToRead = Math.min(molecule.num_atoms - linesAlreadyRead, chunkSize);
        for (let i = 0; i < linesToRead; i++) {
            let subst_idWritten = false;
            let subst_nameWritten = false;
            let chargeWritten = false;
            let status_bitWritten = false;
            for (let j = 0; j < numOfColumn; j++) {
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
                    default:
                        if (hasSubst_id === true && subst_idWritten === false) {
                            TokenBuilder.addUnchecked(subst_idTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                            subst_idWritten = true;
                        } else if (hasSubst_name === true && subst_nameWritten === false) {
                            TokenBuilder.addUnchecked(subst_nameTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                            subst_nameWritten = true;
                        } else if (hasCharge === true && chargeWritten === false) {
                            TokenBuilder.addUnchecked(chargeTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                            chargeWritten = true;
                        } else if (hasStatus_bit === true && status_bitWritten === false) {
                            TokenBuilder.addUnchecked(status_bitTokens, tokenizer.tokenStart, tokenizer.tokenEnd);
                            status_bitWritten = true;
                        }
                }
            }
        }
        linesAlreadyRead += linesToRead;
        return linesToRead;
    }, ctx => ctx.update({ message: 'Parsing...', current: tokenizer.position, max: length }));

    const ret = {
        count: molecule.num_atoms,
        atom_id: atom_idTokenColumn(Column.Schema.int),
        atom_name: atom_nameTokenColumn(Column.Schema.str),
        x: xTokenColumn(Column.Schema.float),
        y: yTokenColumn(Column.Schema.float),
        z: zTokenColumn(Column.Schema.float),
        atom_type: atom_typeColumn(Column.Schema.str),

        // optional columns
        subst_id: hasSubst_id ? subst_idTokenColumn(Column.Schema.int) : undefInt,
        subst_name: hasSubst_name ? subst_nameTokenColumn(Column.Schema.str) : undefStr,
        charge: hasCharge ? chargeTokenColumn(Column.Schema.float) : undefFloat,
        status_bit: hasStatus_bit ? status_bitTokenColumn(Column.Schema.str) : undefStr,
    };
    return ret;
}

async function handleBonds(state: State): Promise<Schema.Mol2Bonds> {
    const { tokenizer, molecule } = state;
    let hasStatus_bit = false;

    while (getTokenString(tokenizer) !== '@<TRIPOS>BOND' && tokenizer.position < tokenizer.data.length) {
        markLine(tokenizer);
    }

    const initialTokenizerPosition = tokenizer.position;
    const initialTokenizerLineNumber = tokenizer.lineNumber;
    const firstLine = readLine(tokenizer);
    const firstLineArray = firstLine.trim().split(/\s+/g);
    const firstLineLength = firstLineArray.length;
    if (firstLineLength === 5) {
        hasStatus_bit = true;
    }

    // required columns
    const bond_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const origin_bond_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const target_bond_idTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const bondTypeTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);

    const bond_idTokenColumn = TokenColumn(bond_idTokens);
    const origin_bond_idTokenColumn = TokenColumn(origin_bond_idTokens);
    const target_bond_idTokenColumn = TokenColumn(target_bond_idTokens);
    const bondTypeTokenColumn = TokenColumn(bondTypeTokens);

    // optional columns
    const status_bitTokens = TokenBuilder.create(tokenizer.data, molecule.num_bonds * 2);
    const status_bitTokenColumn = TokenColumn(status_bitTokens);
    const undefStr = Column.Undefined(molecule.num_bonds, Column.Schema.str);

    let numberOfColumn = 4;
    if (hasStatus_bit) { numberOfColumn++; }

    tokenizer.position = initialTokenizerPosition;
    tokenizer.lineNumber = initialTokenizerLineNumber;

    const { length } = tokenizer;
    let linesAlreadyRead = 0;
    await chunkedSubtask(state.runtimeCtx, 100000, void 0, chunkSize => {
        const linesToRead = Math.min(molecule.num_bonds - linesAlreadyRead, chunkSize);
        for (let i = 0; i < linesToRead; i++) {
            for (let j = 0; j < numberOfColumn; j++) {
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
        bond_id: bond_idTokenColumn(Column.Schema.int),
        origin_atom_id: origin_bond_idTokenColumn(Column.Schema.int),
        target_atom_id: target_bond_idTokenColumn(Column.Schema.int),
        bond_type: bondTypeTokenColumn(Column.Schema.str),

        // optional columns
        status_bits: hasStatus_bit ? status_bitTokenColumn(Column.Schema.str) : undefStr,
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
    }

    const result: Schema.Mol2File = { name, structures };
    return Result.success(result);
}

export function parseMol2(data: string, name: string) {
    return Task.create<Result<Schema.Mol2File>>('Parse MOL2', async ctx => {
        return await parseInternal(ctx, data, name);
    });
}