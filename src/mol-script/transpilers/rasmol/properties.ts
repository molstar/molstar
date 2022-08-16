/*
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
*/

import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import { PropertyDict } from '../types';

const reFloat = /[-+]?[0-9]*\.?[0-9]+/;
const rePosInt = /[0-9]+/;

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

export const properties: PropertyDict = {
    atomno: {
        '@desc': 'sequential number; you can use "@" instead of "atomno=" -- for example, select @33 or Var x = @33 or @35',
        '@examples': ['atomno = 10'],
        isNumeric: true,
        regex: rePosInt, map: x => parseInt(x),
        level: 'atom-test', property: B.ammp('id')
    },
    atomX: {
        '@desc': 'Cartesian X coordinate (or just X)',
        '@examples': ['x = 4.2'],
        abbr: ['X'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('x')
    },
    atomY: {
        '@desc': 'Cartesian Y coordinate (or just Y)',
        '@examples': ['y < 42'],
        abbr: ['Y'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('y')
    },
    atomZ: {
        '@desc': 'Cartesian Z coordinate (or just Z)',
        '@examples': ['Z > 10'],
        abbr: ['Z'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('z')
    },
    elemno: {
        '@desc': 'atomic element number',
        '@examples': ['elemno=8'],
        isNumeric: true,
        regex: /[0-9\s{}-]+/, map: x => parseInt(x),
        level: 'atom-test', property: B.acp('atomicNumber')
    },
    formalCharge: {
        '@desc': 'formal charge',
        '@examples': ['formalCharge=1'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.ammp('pdbx_formal_charge')
    },
    mass: {
        '@desc': 'atomic mass -- especially useful with appended .max or .sum',
        '@examples': ['mass > 13'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('mass')
    },
    occupancy: {
        '@desc': 'CIF file site occupancy. In SELECT command comparisons ("select occupancy < 90"), an integer n implies measurement on a 0-100 scale; also, in the context %[occupancy] or %q for a label, the reported number is a percentage. In all other cases, such as when %Q is used in a label or when a decimal number is used in a comparison, the scale is 0.0 - 1.0.',
        '@examples': ['occupancy < 1'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.ammp('occupancy')
    },
    resno: {
        '@desc': 'PDB residue number, not including insertion code (see also seqcode, below)',
        '@examples': ['resno = 100'],
        isNumeric: true,
        regex: /-?[0-9]+/, map: x => parseInt(x),
        level: 'residue-test', property: B.ammp('auth_seq_id')
    },
    temperature: {
        '@desc': 'yes  yes  temperature factor (B-factor)',
        '@examples': ['temperature >= 20'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.ammp('B_iso_or_equiv')
    },
    vanderwaals: {
        '@desc': 'van der Waals radius',
        '@examples': ['vanderwaals >2'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('vdw')
    },
};
