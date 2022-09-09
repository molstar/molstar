/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * Adapted from MolQL project
 */

import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import { PropertyDict } from '../types';

const reFloat = /[-+]?[0-9]*\.?[0-9]+/;
const rePosInt = /[+]?[0-9]+/;
const reInt = /[-+]?[0-9]+/;

function str(x: string) { return x; }

export const sstrucDict: { [key: string]: string } = {
    T: 'turn', // Turn
    E: 'sheet', // Extended conformation ($\beta$ sheets)
    B: 'strand', // Isolated bridge
    H: 'alpha', // Alpha helix
    G: '3-10', // 3-10 helix
    I: 'pi', // Pi helix
    C: 'none', // Coil
};
export function sstrucMap(x: string) {
    return B.struct.type.secondaryStructureFlags(
        [sstrucDict[x.toUpperCase()] || 'none']
    );
}

export const properties: PropertyDict = {
    name: {
        '@desc': 'str    atom name',
        '@examples': ['name CA'],
        regex: /[a-zA-Z0-9]+/, map: B.atomName,
        level: 'atom-test', property: B.ammp('label_atom_id')
    },
    type: {
        '@desc': 'str    atom type',
        '@examples': ['type C3'],
        isUnsupported: true,
        regex: /[a-zA-Z0-9]+/, map: str,
        level: 'atom-test'
    },
    index: {
        '@desc': 'num    the atom number, starting at 0',
        '@examples': ['index 10'],
        isNumeric: true,
        regex: rePosInt, map: x => (parseInt(x) - 1),
        level: 'atom-test', property: B.ammp('id')
    },
    serial: {
        '@desc': 'num    the atom number, starting at 1',
        '@examples': ['serial 11'],
        isNumeric: true,
        regex: rePosInt, map: x => parseInt(x),
        level: 'atom-test', property: B.ammp('id')
    },
    atomicnumber: {
        '@desc': 'num    atomic number (0 if undefined)',
        '@examples': ['atomicnumber 13'],
        isNumeric: true,
        regex: rePosInt, map: x => parseInt(x),
        level: 'atom-test', property: B.acp('atomicNumber')
    },
    element: {
        '@desc': 'str  atomic element symbol string ("X" if undefined)',
        '@examples': ['element N'],
        regex: /[a-zA-Z0-9]{1,3}/, map: x => B.es(x),
        level: 'atom-test', property: B.acp('elementSymbol')
    },
    altloc: {
        '@desc': 'str  alternate location/conformation identifier',
        '@examples': ['altloc C'],
        regex: /[a-zA-Z0-9]+/, map: str,
        level: 'atom-test', property: B.ammp('label_alt_id')
    },
    chain: {
        '@desc': 'str  the one-character chain identifier',
        '@examples': ['chain A'],
        regex: /[a-zA-Z0-9]+/, map: str,
        level: 'residue-test', property: B.ammp('auth_asym_id')
    },
    residue: {
        '@desc': 'num  a set of connected atoms with the same residue number',
        '@examples': ['residue < 11', 'residue 11'],
        isNumeric: true,
        regex: reInt, map: x => parseInt(x),
        level: 'residue-test', property: B.ammp('auth_seq_id')
    },
    fragment: {
        '@desc': 'num  a set of connected residues',
        '@examples': ['fragment 42'],
        isUnsupported: true,
        isNumeric: true,
        regex: reInt, map: x => parseInt(x),
        level: 'residue-test'
    },
    pfrag: {
        '@desc': 'num  a set of connected protein residues',
        '@examples': ['pfrag 42'],
        isUnsupported: true,
        isNumeric: true,
        regex: reInt, map: x => parseInt(x),
        level: 'residue-test'
    },
    nfrag: {
        '@desc': 'num  a set of connected nucleic residues',
        '@examples': ['nfrag 42'],
        isUnsupported: true,
        isNumeric: true,
        regex: reInt, map: x => parseInt(x),
        level: 'residue-test'
    },
    sequence: {
        '@desc': 'str  a sequence given by one letter names',
        '@examples': ['sequence PGATTACA'],
        isUnsupported: true,
        regex: /[a-zA-Z0-9]+/, map: str,
        level: 'residue-test'
    },
    numbonds: {
        '@desc': 'num  number of bonds',
        '@examples': ['numbonds = 2', 'numbonds >= 3'],
        isNumeric: true,
        regex: rePosInt, map: x => parseInt(x),
        level: 'atom-test', property: B.acp('bondCount')
    },
    resname: {
        '@desc': 'str  residue name',
        '@examples': ['resname ALA'],
        regex: /[a-zA-Z0-9]+/, map: str,
        level: 'residue-test', property: B.ammp('auth_comp_id')
    },
    resid: {
        '@desc': 'num  residue id',
        '@examples': ['resid 42'],
        isNumeric: true,
        regex: reInt, map: x => parseInt(x),
        level: 'residue-test', property: B.ammp('auth_seq_id')
    },
    segname: {
        '@desc': 'str  segment name',
        '@examples': ['segname B'],
        regex: /[a-zA-Z0-9]+/, map: str,
        level: 'residue-test', property: B.ammp('label_asym_id')
    },
    x: {
        '@desc': 'float  x coordinate',
        '@examples': ['x 42'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('x')
    },
    y: {
        '@desc': 'float  y coordinate',
        '@examples': ['y > 1.7'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('y')
    },
    z: {
        '@desc': 'float  z coordinate',
        '@examples': ['z < 11', 'z > -21'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('z')
    },
    radius: {
        '@desc': 'float  atomic radius',
        '@examples': ['radius > 1.3'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('vdw')
    },
    mass: {
        '@desc': 'float  atomic mass',
        '@examples': ['mass > 2'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('mass')
    },
    charge: {
        '@desc': 'float  atomic charge',
        '@examples': ['charge > 0', 'charge 1'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.ammp('pdbx_formal_charge')
    },
    beta: {
        '@desc': 'float  temperature factor',
        '@examples': ['beta < 20', 'beta > 35'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.ammp('B_iso_or_equiv')
    },
    occupancy: {
        '@desc': 'float  occupancy',
        '@examples': ['occupancy 1', 'occupancy < 1'],
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.ammp('occupancy')
    },
    user: {
        '@desc': 'float  time-varying user-specified value',
        '@examples': ['user < 0.1'],
        isUnsupported: true,
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test'
    },
    rasmol: {
        '@desc': 'str  translates Rasmol selection string to VMD',
        '@examples': ["rasmol 'all'"],
        isUnsupported: true,
        regex: /[^']*/, map: str,
        level: 'atom-test'
    },
    structure: {
        '@desc': 'str  single letter name for the secondary structure',
        '@examples': ['structure H', 'structure H E'],
        regex: /T|E|B|H|G|I|C/i, map: sstrucMap,
        level: 'atom-test', property: B.ammp('secondaryStructureFlags')
    },
    phi: {
        '@desc': 'float  phi backbone conformational angles',
        '@examples': ['phi < 160'],
        isUnsupported: true,
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'residue-test'
    },
    psi: {
        '@desc': 'float  psi backbone conformational angles',
        '@examples': ['psi < 160'],
        isUnsupported: true,
        isNumeric: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'residue-test'
    },
    ufx: {
        '@desc': 'num  force to apply in the x coordinate',
        '@examples': ['ufx 1'],
        isUnsupported: true,
        isNumeric: true,
        regex: reFloat, map: x => parseInt(x),
        level: 'atom-test'
    },
    ufy: {
        '@desc': 'num  force to apply in the y coordinate',
        '@examples': ['ufy 1'],
        isUnsupported: true,
        isNumeric: true,
        regex: reFloat, map: x => parseInt(x),
        level: 'atom-test'
    },
    ufz: {
        '@desc': 'num  force to apply in the z coordinate',
        '@examples': ['ufz 1'],
        isUnsupported: true,
        isNumeric: true,
        regex: reFloat, map: x => parseInt(x),
        level: 'atom-test'
    },
};
