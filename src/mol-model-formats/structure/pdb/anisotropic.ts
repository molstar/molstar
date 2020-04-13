/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifField } from '../../../mol-io/reader/cif';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { TokenBuilder, Tokenizer } from '../../../mol-io/reader/common/text/tokenizer';
import { guessElementSymbolTokens } from '../util';
import { parseIntSkipLeadingWhitespace as fastParseInt } from '../../../mol-io/reader/common/text/number-parser';

type AnisotropicTemplate = typeof getAnisotropicTemplate extends (...args: any) => infer T ? T : never
export function getAnisotropicTemplate(data: string, count: number) {
    const str = () => [] as string[];
    const float = () => new Float32Array(count);
    const ts = () => TokenBuilder.create(data, 2 * count);
    return {
        index: 0,
        count,
        id: str(),
        type_symbol: ts(),
        pdbx_label_atom_id: ts(),
        pdbx_label_alt_id: ts(),
        pdbx_label_comp_id: ts(),
        pdbx_label_asym_id: ts(),
        pdbx_label_seq_id: ts(),
        pdbx_PDB_ins_code: ts(),
        'U[1][1]': float(),
        'U[2][2]': float(),
        'U[3][3]': float(),
        'U[1][2]': float(),
        'U[1][3]': float(),
        'U[2][3]': float(),
        pdbx_auth_seq_id: ts(),
        pdbx_auth_comp_id: ts(),
        pdbx_auth_asym_id: ts(),
        pdbx_auth_atom_id: ts(),
    };
}

export function getAnisotropic(sites: AnisotropicTemplate): { [K in keyof mmCIF_Schema['atom_site_anisotrop']]?: CifField } {
    const pdbx_auth_seq_id = CifField.ofTokens(sites.pdbx_auth_seq_id);
    const pdbx_auth_comp_id = CifField.ofTokens(sites.pdbx_auth_comp_id);
    const pdbx_auth_asym_id = CifField.ofTokens(sites.pdbx_auth_asym_id);
    const pdbx_auth_atom_id = CifField.ofTokens(sites.pdbx_auth_atom_id);

    const fields: { [K in keyof mmCIF_Schema['atom_site_anisotrop']]?: CifField } = {
        id: CifField.ofStrings(sites.id),
        type_symbol: CifField.ofTokens(sites.type_symbol),
        pdbx_label_atom_id: pdbx_auth_atom_id,
        pdbx_label_alt_id: CifField.ofTokens(sites.pdbx_label_alt_id),
        pdbx_label_comp_id: pdbx_auth_comp_id,
        pdbx_label_asym_id: pdbx_auth_asym_id,
        pdbx_label_seq_id: pdbx_auth_seq_id,
        pdbx_PDB_ins_code: CifField.ofTokens(sites.pdbx_PDB_ins_code),

        pdbx_auth_seq_id,
        pdbx_auth_comp_id,
        pdbx_auth_asym_id,
        pdbx_auth_atom_id,
    };

    (fields as any)['U[1][1]'] = CifField.ofNumbers(sites['U[1][1]']);
    (fields as any)['U[2][2]'] = CifField.ofNumbers(sites['U[2][2]']);
    (fields as any)['U[3][3]'] = CifField.ofNumbers(sites['U[3][3]']);
    (fields as any)['U[1][2]'] = CifField.ofNumbers(sites['U[1][2]']);
    (fields as any)['U[1][3]'] = CifField.ofNumbers(sites['U[1][3]']);
    (fields as any)['U[2][3]'] = CifField.ofNumbers(sites['U[2][3]']);

    return fields;
}

export function addAnisotropic(sites: AnisotropicTemplate, model: string, data: Tokenizer, s: number, e: number) {
    const { data: str } = data;
    const length = e - s;

    // COLUMNS       DATA  TYPE    FIELD          DEFINITION
    // -----------------------------------------------------------------
    // 1 - 6        Record name   "ANISOU"
    // 7 - 11       Integer       serial         Atom serial number.
    Tokenizer.trim(data, s + 6, s + 11);
    sites.id[sites.index] = str.substring(data.tokenStart, data.tokenEnd);

    // 13 - 16       Atom          name           Atom name.
    TokenBuilder.addToken(sites.pdbx_auth_atom_id, Tokenizer.trim(data, s + 12, s + 16));

    // 17            Character     altLoc         Alternate location indicator
    if (str.charCodeAt(s + 16) === 32) { // ' '
        TokenBuilder.add(sites.pdbx_label_alt_id, 0, 0);
    } else {
        TokenBuilder.add(sites.pdbx_label_alt_id, s + 16, s + 17);
    }

    // 18 - 20       Residue name  resName        Residue name.
    TokenBuilder.addToken(sites.pdbx_auth_comp_id, Tokenizer.trim(data, s + 17, s + 20));

    // 22            Character     chainID        Chain identifier.
    TokenBuilder.add(sites.pdbx_auth_asym_id, s + 21, s + 22);

    // 23 - 26       Integer       resSeq         Residue sequence number.
    TokenBuilder.addToken(sites.pdbx_auth_seq_id, Tokenizer.trim(data, s + 22, s + 26));

    // 27            AChar         iCode          Insertion code.
    if (str.charCodeAt(s + 26) === 32) { // ' '
        TokenBuilder.add(sites.pdbx_PDB_ins_code, 0, 0);
    } else {
        TokenBuilder.add(sites.pdbx_PDB_ins_code, s + 26, s + 27);
    }

    // 29 - 35       Integer       u[0][0]        U(1,1)
    sites['U[1][1]'][sites.index] = fastParseInt(str, s + 28, s + 35) / 10000;

    // 36 - 42       Integer       u[1][1]        U(2,2)
    sites['U[2][2]'][sites.index] = fastParseInt(str, s + 35, s + 42) / 10000;

    // 43 - 49       Integer       u[2][2]        U(3,3)
    sites['U[3][3]'][sites.index] = fastParseInt(str, s + 42, s + 49) / 10000;

    // 50 - 56       Integer       u[0][1]        U(1,2)
    sites['U[1][2]'][sites.index] = fastParseInt(str, s + 49, s + 56) / 10000;

    // 57 - 63       Integer       u[0][2]        U(1,3)
    sites['U[1][3]'][sites.index] = fastParseInt(str, s + 56, s + 63) / 10000;

    // 64 - 70       Integer       u[1][2]        U(2,3)
    sites['U[2][3]'][sites.index] = fastParseInt(str, s + 63, s + 70) / 10000;

    // 77 - 78       LString(2)    element        Element symbol, right-justified.
    if (length >= 78) {
        Tokenizer.trim(data, s + 76, s + 78);

        if (data.tokenStart < data.tokenEnd) {
            TokenBuilder.addToken(sites.type_symbol, data);
        } else {
            guessElementSymbolTokens(sites.type_symbol, str, s + 12, s + 16);
        }
    } else {
        guessElementSymbolTokens(sites.type_symbol, str, s + 12, s + 16);
    }

    // 79 - 80       LString(2)    charge         Charge on the atom.
    // TODO

    sites.index++;
}