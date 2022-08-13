/*
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * @author Koya Sakuma
 * This module was taken from jmol transpiler from MolQL and modified in similar manner as pymol and vmd tranpilers.                                             \
*/

import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import { PropertyDict } from '../types';

const reFloat = /[-+]?[0-9]*\.?[0-9]+/;
const rePosInt = /[0-9]+/;

function str(x: string) { return x; }

export function sstrucMap(x: string) {
    return B.struct.type.secondaryStructureFlags(
        [structureDict[x.toUpperCase()] || 'none']
    );
}


export const structureDict: {[key: string]: string} = {
    none: 'none',
    turn: 'turn',
    sheet: 'beta',
    helix: 'helix',
    dna: 'dna',
    rna: 'rna',

    0: 'none',
    1: 'turn',
    2: 'beta',
    3: 'helix',
    4: 'dna',
    5: 'rna',
    6: 'carbohydrate',
    7: '3-10',
    8: 'alpha',
    9: 'pi',
};
export function structureMap(x: any) {
    if (x.head && x.head === 'core.type.regex') x = x.args[0].replace(/^\^|\$$/g, '');
    x = structureDict[x.toString().toLowerCase()] || 'none';
    if (['dna', 'rna', 'carbohydrate'].indexOf(x) !== -1) {
        throw new Error("values 'dna', 'rna', 'carbohydrate' not yet supported for 'structure' property");
    } else {
        return B.struct.type.secondaryStructureFlags([x]);
    }
}

export const special: PropertyDict = {
    hoge: {
        '@desc': 'PDB residue number, not including insertion code (see also seqcode, below)',
        '@examples': ['resno = 100'],
//        isNumeric: true,
        regex: /-?[0-9]+/, map: x => parseInt(x),
        level: 'residue-test', property: B.ammp('auth_seq_id')
    },
};

