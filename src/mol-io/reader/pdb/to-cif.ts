/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifField, CifCategory } from '../cif';
import { mmCIF_Schema } from '../cif/schema/mmcif';
import CifTextField, { CifTextValueField } from '../cif/text/field';
import { TokenBuilder, Tokenizer } from '../common/text/tokenizer';
import { PdbFile } from './schema';
import { CifFile } from '../cif/data-model';
import { substringStartsWith } from 'mol-util/string';
import { Task } from 'mol-task';

function toCategory(name: string, fields: { [name: string]: CifField | undefined }, rowCount: number): CifCategory {
    return {
        name,
        fieldNames: Object.keys(fields),
        rowCount,
        getField(f: string) {
            return fields[f];
        }
    }
}

function _entity(): { [K in keyof mmCIF_Schema['entity']]?: CifField } {
    return {
        id: CifTextValueField(['1', '2', '3']),
        type: CifTextValueField(['polymer', 'non-polymer', 'water'])
    }
}

function atom_site_template(data: string, count: number) {
    const str = () => new Array(count) as string[];
    const ts = () => TokenBuilder.create(data, 2 * count);
    return {
        index: 0,
        count,
        group_PDB: ts(),
        id: str(),
        auth_atom_id: ts(),
        label_alt_id: ts(),
        auth_comp_id: ts(),
        auth_asym_id: ts(),
        auth_seq_id: ts(),
        pdbx_PDB_ins_code: ts(),
        Cartn_x: ts(),
        Cartn_y: ts(),
        Cartn_z: ts(),
        occupancy: ts(),
        B_iso_or_equiv: ts(),
        type_symbol: ts(),
        pdbx_PDB_model_num: str(),
        label_entity_id: str()
    };
}

function _atom_site(sites: AtomSiteTemplate): { [K in keyof mmCIF_Schema['atom_site']]?: CifField } {
    const auth_asym_id = CifTextField(sites.auth_asym_id, sites.count);
    const auth_atom_id = CifTextField(sites.auth_atom_id, sites.count);
    const auth_comp_id = CifTextField(sites.auth_comp_id, sites.count);
    const auth_seq_id = CifTextField(sites.auth_seq_id, sites.count);

    return {
        auth_asym_id,
        auth_atom_id,
        auth_comp_id,
        auth_seq_id,
        B_iso_or_equiv: CifTextField(sites.B_iso_or_equiv, sites.count),
        Cartn_x: CifTextField(sites.Cartn_x, sites.count),
        Cartn_y: CifTextField(sites.Cartn_y, sites.count),
        Cartn_z: CifTextField(sites.Cartn_z, sites.count),
        group_PDB: CifTextField(sites.group_PDB, sites.count),
        id: CifTextValueField(sites.id),

        label_alt_id: CifTextField(sites.label_alt_id, sites.count),

        label_asym_id: auth_asym_id,
        label_atom_id: auth_atom_id,
        label_comp_id: auth_comp_id,
        label_seq_id: auth_seq_id,
        label_entity_id: CifTextValueField(sites.label_entity_id),

        occupancy: CifTextField(sites.occupancy, sites.count),
        type_symbol: CifTextField(sites.type_symbol, sites.count),

        pdbx_PDB_ins_code: CifTextField(sites.pdbx_PDB_ins_code, sites.count),
        pdbx_PDB_model_num: CifTextValueField(sites.pdbx_PDB_model_num)
    };
}

const WaterNames = new Set([ 'SOL', 'WAT', 'HOH', 'H2O', 'W', 'DOD', 'D3O', 'TIP3', 'TIP4', 'SPC' ]);

function getEntityId(residueName: string, isHet: boolean) {
    if (isHet) {
        if (WaterNames.has(residueName)) return '3';
        return '2';
    }
    return '1';
}

function addAtom(sites: AtomSiteTemplate, model: string, data: Tokenizer, s: number, e: number, isHet: boolean) {
    const { data: str } = data;
    let startPos = s;
    let start = s;
    const end = e;
    const length = end - start;

    // TODO: filter invalid atoms

    // COLUMNS        DATA TYPE       CONTENTS
    // --------------------------------------------------------------------------------
    // 1 -  6        Record name     "ATOM  "
    Tokenizer.trim(data, start, start + 6);
    TokenBuilder.add(sites.group_PDB, data.tokenStart, data.tokenEnd);

    // 7 - 11        Integer         Atom serial number.
    // TODO: support HEX
    start = startPos + 6;
    Tokenizer.trim(data, start, start + 5);
    sites.id[sites.index] = data.data.substring(data.tokenStart, data.tokenEnd);

    // 13 - 16        Atom            Atom name.
    start = startPos + 12;
    Tokenizer.trim(data, start, start + 4);
    TokenBuilder.add(sites.auth_atom_id, data.tokenStart, data.tokenEnd);

    // 17             Character       Alternate location indicator.
    if (str.charCodeAt(startPos + 16) === 32) { // ' '
        TokenBuilder.add(sites.label_alt_id, 0, 0);
    } else {
        TokenBuilder.add(sites.label_alt_id, startPos + 16, startPos + 17);
    }

    // 18 - 20        Residue name    Residue name.
    start = startPos + 17;
    Tokenizer.trim(data, start, start + 3);
    TokenBuilder.add(sites.auth_comp_id, data.tokenStart, data.tokenEnd);
    const residueName = str.substring(data.tokenStart, data.tokenEnd);

    // 22             Character       Chain identifier.
    TokenBuilder.add(sites.auth_asym_id, startPos + 21, startPos + 22);

    // 23 - 26        Integer         Residue sequence number.
    // TODO: support HEX
    start = startPos + 22;
    Tokenizer.trim(data, start, start + 4);
    TokenBuilder.add(sites.auth_seq_id, data.tokenStart, data.tokenEnd);

    // 27             AChar           Code for insertion of residues.
    if (str.charCodeAt(startPos + 26) === 32) { // ' '
        TokenBuilder.add(sites.label_alt_id, 0, 0);
    } else {
        TokenBuilder.add(sites.label_alt_id, startPos + 26, startPos + 27);
    }

    // 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
    start = startPos + 30;
    Tokenizer.trim(data, start, start + 8);
    TokenBuilder.add(sites.Cartn_x, data.tokenStart, data.tokenEnd);

    // 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
    start = startPos + 38;
    Tokenizer.trim(data, start, start + 8);
    TokenBuilder.add(sites.Cartn_y, data.tokenStart, data.tokenEnd);

    // 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
    start = startPos + 46;
    Tokenizer.trim(data, start, start + 8);
    TokenBuilder.add(sites.Cartn_z, data.tokenStart, data.tokenEnd);

    // 55 - 60        Real(6.2)       Occupancy.
    start = startPos + 54;
    Tokenizer.trim(data, start, start + 6);
    TokenBuilder.add(sites.occupancy, data.tokenStart, data.tokenEnd);

    // 61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
    if (length >= 66) {
        start = startPos + 60;
        Tokenizer.trim(data, start, start + 6);
        TokenBuilder.add(sites.B_iso_or_equiv, data.tokenStart, data.tokenEnd);
    } else {
        TokenBuilder.add(sites.label_alt_id, 0, 0);
    }

    // 73 - 76        LString(4)      Segment identifier, left-justified.
    // ignored

    // 77 - 78        LString(2)      Element symbol, right-justified.
    if (length >= 78) {
        start = startPos + 76;
        Tokenizer.trim(data, start, start + 2);

        if (data.tokenStart < data.tokenEnd) {
            TokenBuilder.add(sites.type_symbol, data.tokenStart, data.tokenEnd);
        } else {
            // "guess" the symbol
            TokenBuilder.add(sites.type_symbol, startPos + 12, startPos + 13);
        }
    } else {
        TokenBuilder.add(sites.type_symbol, startPos + 12, startPos + 13);
    }

    sites.label_entity_id[sites.index] = getEntityId(residueName, isHet);
    sites.pdbx_PDB_model_num[sites.index] = model;

    sites.index++;
}

type AtomSiteTemplate = typeof atom_site_template extends (...args: any) => infer T ? T : never

async function pdbToMmCIF(pdb: PdbFile): Promise<CifFile> {
    const { lines } = pdb;
    const { data, indices } = lines;
    const tokenizer = Tokenizer(data);

    // Count the atoms
    let atomCount = 0;
    for (let i = 0, _i = lines.count; i < _i; i++) {
        const s = indices[2 * i], e = indices[2 * i + 1];
        switch (data[s]) {
            case 'A':
                if (substringStartsWith(data, s, e, 'ATOM  ')) atomCount++;
                break;
            case 'H':
                if (substringStartsWith(data, s, e, 'HETATM')) atomCount++;
                break;
        }
    }

    const atom_site = atom_site_template(data, atomCount);

    let modelNum = 0, modelStr = '';

    for (let i = 0, _i = lines.count; i < _i; i++) {
        const s = indices[2 * i], e = indices[2 * i + 1];
        switch (data[s]) {
            case 'A':
                if (!substringStartsWith(data, s, e, 'ATOM  ')) continue;
                if (!modelNum) { modelNum++; modelStr = '' + modelNum; }
                addAtom(atom_site, modelStr, tokenizer, s, e, false);
                break;
            case 'H':
                if (!substringStartsWith(data, s, e, 'HETATM')) continue;
                if (!modelNum) { modelNum++; modelStr = '' + modelNum; }
                addAtom(atom_site, modelStr, tokenizer, s, e, true);
                break;
            case 'M':
                if (substringStartsWith(data, s, e, 'MODEL ')) {
                    modelNum++;
                    modelStr = '' + modelNum;
                }
                break;

        }
    }

    const categories = {
        entity: toCategory('entity', _entity(), 3),
        atom_site: toCategory('atom_site', _atom_site(atom_site), atomCount)
    }

    return {
        name: pdb.id,
        blocks: [{
            saveFrames: [],
            header: pdb.id || 'PDB',
            categoryNames: Object.keys(categories),
            categories
        }]
    };
}

export function convertPDBtoMmCif(pdb: PdbFile): Task<CifFile> {
    return Task.create('Convert PDB to mmCIF', async ctx => {
        await ctx.update('Converting to mmCIF...');
        return pdbToMmCIF(pdb);
    });
}