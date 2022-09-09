/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Koya Sakuma <koya.sakuma.work@gmail.com>
 *
 * Adapted from MolQL project
 */

import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import { PropertyDict } from '../types';

const reFloat = /[-+]?[0-9]*\.?[0-9]+/;
const rePosInt = /[0-9]+/;

function str(x: string) { return x; }

const structureDict: { [key: string]: string } = {
    none: 'none',
    turn: 'turn',
    sheet: 'beta',
    helix: 'helix',
    dna: 'dna',
    rna: 'rna',
    carbohydrate: 'carbohydrate',
    helix310: '3-10',
    helixalpha: 'alpha',
    helixpi: 'pi',

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
    if (x.head) {
        if (x.head.name && x.head.name === 'core.type.regex') x = x.args[0].replace(/^\^|\$$/g, '');
        x = structureDict[x.toString().toLowerCase()] || 'none';
        if (['dna', 'rna', 'carbohydrate'].indexOf(x) !== -1) {
            throw new Error("values 'dna', 'rna', 'carbohydrate' not yet supported for 'structure' property");
        } else {
            return B.struct.type.secondaryStructureFlags([x]);
        }
    }
}

export const properties: PropertyDict = {
    adpmax: {
        '@desc': 'the maximum anisotropic displacement parameter for the selected atom',
        '@examples': [''],
        isUnsupported: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test'
    },
    adpmin: {
        '@desc': 'the minimum anisotropic displacement parameter for the selected atom',
        '@examples': [''],
        isUnsupported: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test'
    },
    altloc: {
        '@desc': 'PDB alternate location identifier',
        '@examples': ['altloc = A'],
        regex: /[a-zA-Z0-9]/, map: str,
        level: 'atom-test', property: B.ammp('label_alt_id')
    },
    altname: {
        '@desc': 'an alternative name given to atoms by some file readers (for example, P2N)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[a-zA-Z0-9]/, map: str,
        level: 'atom-test'
    },
    atomID: {
        '@desc': 'special atom IDs for PDB atoms assigned by Jmol',
        '@examples': [''],
        isUnsupported: true,
        regex: rePosInt, map: x => parseInt(x),
        level: 'atom-test'
    },
    atomIndex: {
        '@desc': 'atom 0-based index; a unique number for each atom regardless of the number of models loaded',
        '@examples': [''],
        isUnsupported: true,
        regex: rePosInt, map: x => parseInt(x),
        level: 'atom-test'
    },
    atomName: {
        '@desc': 'atom name',
        '@examples': ['atomName = CA'],
        regex: /[a-zA-Z0-9]+/, map: v => B.atomName(v),
        level: 'atom-test', property: B.ammp('label_atom_id')
    },
    atomno: {
        '@desc': 'sequential number; you can use "@" instead of "atomno=" -- for example, select @33 or Var x = @33 or @35',
        '@examples': [''],
        isUnsupported: true,
        regex: rePosInt, map: x => parseInt(x),
        level: 'atom-test'
    },
    atomType: {
        '@desc': 'atom type (mol2, AMBER files) or atom name (other file types)',
        '@examples': ['atomType = OH'],
        regex: /[a-zA-Z0-9]+/, map: v => B.atomName(v),
        level: 'atom-test', property: B.ammp('label_atom_id')
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
    bondcount: {
        '@desc': 'covalent bond count',
        '@examples': ['bondcount = 0'],
        isNumeric: true,
        regex: rePosInt, map: x => parseInt(x),
        level: 'atom-test', property: B.acp('bondCount')
    },
    bondingRadius: {
        '@desc': 'radius used for auto bonding; synonymous with ionic and ionicRadius',
        '@examples': [''],
        abbr: ['ionic', 'ionicRadius'],
        isUnsupported: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test'
    },
    cell: {
        '@desc': 'crystallographic unit cell, expressed either in lattice integer notation (111-999) or as a coordinate in ijk space, where {1 1 1} is the same as 555. ANDing two cells, for example select cell=555 and cell=556, selects the atoms on the common face. (Note: in the specifc case of CELL, only "=" is allowed as a comparator.)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    configuration: {
        '@desc': 'Only in the context {configuration=n}, this option selects the set of atoms with either no ALTLOC specified or those atoms having this index into the array of altlocs within its model. So, for example, if the model has altloc "A" and "B", select configuration=1 is equivalent to select altloc="" or altloc="A", and print {configuration=2} is equivalent to print {altloc="" or altloc="B"}. Configuration 0 is "all atoms in a model having configurations", and an invalid configuration number gives no atoms. (Note: in the specifc case of CONFIGURATION, only "=" is allowed as a comparator.)',
        '@examples': [''],
        isUnsupported: true,
        regex: rePosInt, map: x => parseInt(x),
        level: 'atom-test'
    },
    chain: {
        '@desc': 'protein chain. For newer CIF files allowing multicharacter chain specifications, use quotations marks: select chain="AA". For these multicharacter desigations, case is not checked unless the CIF file has lower-case chain designations.',
        '@examples': ['chain = A', 'chain = "AA"'],
        regex: /[a-zA-Z0-9]+/, map: str,
        level: 'chain-test', property: B.ammp('auth_asym_id')
    },
    chainNo: {
        '@desc': 'chain number; sequentially counted from 1 for each model; chainNo == 0 means"no chain" or PDB chain identifier indicated as a blank (Jmol 14.0).',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    color: {
        '@desc': 'the atom color',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    covalentRadius: {
        '@desc': 'covalent bonding radius, synonymous with covalent. Not used by Jmol, but could be used, for example, in {*}.spacefill={*}.covalentRadius.all.',
        '@examples': [''],
        abbr: ['covalent'],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    cs: {
        '@desc': 'chemical shift calculated using computational results that include magnetic shielding tensors.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    element: {
        '@desc': 'element symbol. The value of this parameter depends upon the context. Used with select structure=x, x can be either the quoted element symbol, "H", "He", "Li", etc. or atomic number. In all other contexts, the value is the element symbol. When the atom is a specific isotope, the string will contain the isotope number -- "13C", for example.',
        '@examples': ['element=Fe'],
        regex: /[a-zA-Z]+/, map: x => B.es(x),
        level: 'atom-test', property: B.acp('elementSymbol')
    },
    elemno: {
        '@desc': 'atomic element number',
        '@examples': ['elemno=8'],
        regex: /[0-9\s{}-]+/, map: x => parseInt(x),
        level: 'atom-test', property: B.acp('atomicNumber')
    },
    eta: {
        '@desc': 'Based on Carlos M. Duarte, Leven M. Wadley, and Anna Marie Pyle, RNA structure comparison, motif search and discovery using a reduced representation of RNA conformational space, Nucleic Acids Research, 2003, Vol. 31, No. 16 4755-4761. The parameter eta is the C4\'[i-1]-P[i]-C4\'[i]-P[i+1] dihedral angle; theta is the P[i]-C4\'[i]-P[i+1]-C4\'[i+1] dihedral angle. Both are measured on a 0-360 degree scale because they are commonly near 180 degrees. Using the commands plot PROPERTIES eta theta resno; select visible;wireframe only one can create these authors\' "RNA worm" graph.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    theta: {
        '@desc': 'Based on Carlos M. Duarte, Leven M. Wadley, and Anna Marie Pyle, RNA structure comparison, motif search and discovery using a reduced representation of RNA conformational space, Nucleic Acids Research, 2003, Vol. 31, No. 16 4755-4761. The parameter eta is the C4\'[i-1]-P[i]-C4\'[i]-P[i+1] dihedral angle; theta is the P[i]-C4\'[i]-P[i+1]-C4\'[i+1] dihedral angle. Both are measured on a 0-360 degree scale because they are commonly near 180 degrees. Using the commands plot PROPERTIES eta theta resno; select visible;wireframe only one can create these authors\' "RNA worm" graph.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    file: {
        '@desc': 'file number containing this atom',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    formalCharge: {
        '@desc': 'formal charge',
        '@examples': ['formalCharge=1'],
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.ammp('pdbx_formal_charge')
    },
    format: {
        '@desc': 'format (label) of the atom.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    fXyz: {
        '@desc': 'fractional XYZ coordinates',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    fX: {
        '@desc': 'fractional X coordinate',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    fY: {
        '@desc': 'fractional Y coordinate',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    fZ: {
        '@desc': 'fractional Z coordinate',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    fuxyz: {
        '@desc': 'fractional XYZ coordinates in the unitcell coordinate system',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    fux: {
        '@desc': 'fractional X coordinate in the unitcell coordinate system',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    fuy: {
        '@desc': 'fractional Y coordinate in the unitcell coordinate system',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    fuz: {
        '@desc': 'fractional Z coordinate in the unit cell coordinate system',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    group: {
        '@desc': '3-letter residue code',
        '@examples': ['group = ALA'],
        regex: /[a-zA-Z0-9]{1,3}/, map: str,
        level: 'residue-test', property: B.ammp('label_comp_id')
    },
    group1: {
        '@desc': 'single-letter residue code (amino acids only)',
        '@examples': ['group1 = G'],
        regex: /[a-zA-Z]/, map: str,
        level: 'residue-test', property: B.ammp('label_comp_id')
    },
    groupID: {
        '@desc': 'group ID number: A unique ID for each amino acid or nucleic acid residue in a PDB file. 0  noGroup 1-5  ALA, ARG, ASN, ASP, CYS 6-10  GLN, GLU, GLY, HIS, ILE 11-15  LEU, LYS, MET, PHE, PRO 16-20  SER, THR, TRP, TYR, VAL 21-23  ASX, GLX, UNK 24-29  A, +A, G, +G, I, +I 30-35  C, +C, T, +T, U, +U Additional unique numbers are assigned arbitrarily by Jmol and cannot be used reproducibly.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    groupindex: {
        '@desc': 'overall group index',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    hydrophobicity: {
        '@desc': 'Aminoacid residue scale of hydrophobicity based on Rose, G. D., Geselowitz, A. R., Lesser, G. J., Lee, R. H., and Zehfus, M. H. (1985). Hydrophobicity of amino acid residues in globular proteins, Science, 229(4716):834-838.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    identify: {
        '@desc': 'for a PDB/mmCIF file, a label such as [ILE]7^1:A.CD1%A/3 #47, which includes the group ([ILE]), residue number with optional insertion code (7^1), chain (:A), atom name (CD1), alternate location if present (%A), PDB model number (/3, for NMR models when one file is loaded; /file.model such as /2.3 if more than one file is loaded), and atom number (#47). For non-PDB data, the information is shorter -- for example, H15/2.1 #6, indicating atom name (H15), full file.model number (/2.1), and atom number (#6). If only a single model is loaded, %[identify] does not include the model number.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    insertion: {
        '@desc': 'protein residue insertion code',
        '@examples': ['insertion=A'],
        regex: /[a-zA-Z0-9]/, map: str,
        level: 'atom-test', property: B.ammp('pdbx_PDB_ins_code')
    },
    label: {
        '@desc': 'current atom label (same as format)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    mass: {
        '@desc': 'atomic mass -- especially useful with appended .max or .sum',
        '@examples': ['mass > 13'],
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('mass')
    },
    model: {
        '@desc': 'model number',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    modelindex: {
        '@desc': 'a unique number for each model, starting with 0 and spanning all models in all files',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    modO: {
        '@desc': 'currently calculated occupancy from modulation (0 to 100; NaN if atom has no occupancy modulation)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    modXYZ: {
        '@desc': 'currently calculated displacement modulation (for incommensurately modulated structures). Also modX, modY, modZ for individual components. For atoms without modultion, {xx}.modXYZ is -1 and {xx}.modX is NaN, and in a label %[modXYZ] and %[modX] are blank.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    molecule: {
        '@desc': 'molecule number',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    monomer: {
        '@desc': 'monomer number (group number) in a polymer (usually a chain), starting with 1, or 0 if not part of a biopolymer -- that is, not a connected carbohydrate, amino acid, or nucleic acid (Jmol 14.3.15)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    ms: {
        '@desc': 'magnetic shielding calculated from file-loaded tensors.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    occupancy: {
        '@desc': 'CIF file site occupancy. In SELECT command comparisons ("select occupancy < 90"), an integer n implies measurement on a 0-100 scale; also, in the context %[occupancy] or %q for a label, the reported number is a percentage. In all other cases, such as when %Q is used in a label or when a decimal number is used in a comparison, the scale is 0.0 - 1.0.',
        '@examples': ['occupancy < 1'],
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.ammp('occupancy')
    },
    partialCharge: {
        '@desc': 'partial charge',
        '@examples': [''],
        isUnsupported: true,
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test'
    },
    phi: {
        '@desc': 'protein group PHI angle for atom\'s residue',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    polymer: {
        '@desc': 'sequential polymer number in a model, starting with 1.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    polymerLength: {
        '@desc': 'polymer length',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    property_xx: {
        '@desc': 'a property created using the DATA command',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    psi: {
        '@desc': 'protein group PSI angle for the atom\'s residue',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    radius: {
        '@desc': 'currently displayed radius -- In SELECT command comparisons ("select radius=n"), integer n implies Rasmol units 1/250 Angstroms; in all other cases or when a decimal number is used, the units are Angstroms.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    resno: {
        '@desc': 'PDB residue number, not including insertion code (see also seqcode, below)',
        '@examples': ['resno = 100'],
        regex: /-?[0-9]+/, map: x => parseInt(x),
        level: 'residue-test', property: B.ammp('auth_seq_id')
    },
    selected: {
        '@desc': '1.0 if atom is selected; 0.0 if not',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    sequence: {
        '@desc': 'PDB one-character sequence code, as a string of characters, with "?" indicated where single-character codes are not available',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    seqcode: {
        '@desc': 'PDB residue number, including insertion code (for example, 234^2; "seqcode" option added in Jmol 14.3.16)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    seqid: {
        '@desc': '(mmCIF only) the value from _atom_site.label_seq_id; a pointer to _entity_poly_seq.num in the ENTITY_POLY_SEQ category specifying the sequence of monomers in a polymer. Allowance is made for the possibility of microheterogeneity in a sample by allowing a given sequence number to be correlated with more than one monomer id. (Jmol 14.2.3)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    shape: {
        '@desc': 'hybridization geometry such as "tetrahedral"',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    site: {
        '@desc': 'crystallographic site number',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    spacefill: {
        '@desc': 'currently displayed radius',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    straightness: {
        '@desc': 'quaternion-derived straightness (second derivative of the quaternion describing the orientation of the residue. This quantity will have different values depending upon the setting of quaternionFrame as "A" (alpha-carbon/phosphorus atom only), "C" (alpha-carbon/pyrimidine or purine base based), "P" (carbonyl-carbon peptide plane/phosphorus tetrahedron based), or "N" (amide-nitrogen based). The default is alpha-carbon based, which corresponds closely to the following combination of Ramachandran angles involving three consecutive residues i-1, i, and i+1: -psii-1 - phii + psii + phii+1.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    strucno: {
        '@desc': 'a unique number for each helix, sheet, or turn in a model, starting with 1.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    structure: {
        '@desc': 'The value of this parameter depends upon the context. Used with select structure=x, x can be either the quoted keyword "none", "turn", "sheet", "helix", "dna", "rna", or "carbohydrate" or a respective number 0-6. In the context {*}.structure, the return value is a number; in the context label %[structure], the return is one of the six keywords.',
        '@examples': ['structure="helix"', 'structure=3'],
        regex: /none|turn|sheet|helix|dna|rna|carbohydrate|[0-6]/i, map: str,
        level: 'residue-test', property: 'structure'
    },
    substructure: {
        '@desc': 'like structure, the value of this parameter depends upon the context. Used with select substructure=x, x can be either the quoted keyword "none", "turn", "sheet", "helix", "dna", "rna", "carbohydrate", "helix310", "helixalpha", or "helixpi", or the respective number 0-9. In the context {*}.substructure, the return value is a number; in the context label %[substructure], the return is one of the nine keywords.',
        '@examples': ['substructure = "alphahelix"', 'substructure =9'],
        regex: /none|turn|sheet|helix|dna|rna|carbohydrate|helix310|helixalpha|helixpi|[0-9]/i, map: str,
        level: 'residue-test', property: 'structure'
    },
    surfacedistance: {
        '@desc': 'A value related to the distance of an atom to a nominal molecular surface. 0 indicates at the surface. Positive numbers are minimum distances in Angstroms from the given atom to the surface.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    symop: {
        '@desc': 'the first symmetry operation code that generated this atom by Jmol; an integer starting with 1. See also symmetry, below. This operator is only present if the file contains space group information and the file was loaded using the {i, j, k} option so as to generate symmetry-based atoms. To select only the original atoms prior to application of symmetry, you can either use "SYMOP=n", where n is the symmetry operator corresponding to "x,y,z", or you can specify instead simply "NOT symmetry" the way you might specify "NOT hydrogen". Note that atoms in special positions will have multiple operator matches. These atoms can be selected using the keyword SPECIALPOSITION. The special form select SYMOP=nijk selects a specific translation of atoms from the given crystallographic symmetry operation. Comparators <, <=, >, >=, and != can be used and only apply to the ijk part of the designation. The ijk are relative, not absolute. Thus, symop=2555 selects for atoms that have been transformed by symop=2 but not subjected to any further translation. select symop=1555 is identical to select not symmetry. All other ijk are relative to these selections for 555. If the model was loaded using load "filename.cif" {444 666 1}, where the 1 indicates that all symmetry-generated atoms are to be packed within cell 555 and then translated to fill the other 26 specified cells, then select symop=3555 is nearly the same as select symop=3 and cell=555. (The difference being that cell=555 selects for all atoms that are on any edge of the cell, while symop=3555 does not.) However, the situation is different if instead the model was loaded using load "filename.cif" {444 666 0}, where the 0 indicates that symmetry-generated atoms are to be placed exactly where their symmetry operator would put them (x,-y,z being different then from x, 1-y, z). In that case, select symop=3555 is for all atoms that have been generated using symmetry operation 3 but have not had any additional translations applied to the x,y,z expression found in the CIF file. If, for example, symmetry operation 3 is -x,-y,-z, then load "filename.cif" {444 666 0} will place an atom originally at {1/2, 1/2, 1/2} at positions {-1/2, -1/2, -1/2} (symop=3555) and {-3/2, -3/2, -3/2} (symop=3444) and 24 other sites.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    symmetry: {
        '@desc': 'as "symmetry" or in a label as lower-case "o" gives list of crystallographic symmetry operators generating this atom with lattice designations,such as 3555; upper-case "%O" in a label gives a list without the lattice designations. See also symop, above.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    temperature: {
        '@desc': 'yes  yes  temperature factor (B-factor)',
        '@examples': ['temperature >= 20'],
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.ammp('B_iso_or_equiv')
    },
    unitXyz: {
        '@desc': 'unit cell XYZ coordinates',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    uX: {
        '@desc': 'unit cell X coordinate normalized to [0,1)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    uY: {
        '@desc': 'unit cell Y coordinate normalized to [0,1)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    uZ: {
        '@desc': 'unit cell Z coordinate normalized to [0,1)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    valence: {
        '@desc': 'the valence of an atom (sum of bonds, where double bond counts as 2 and triple bond counts as 3',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    vanderwaals: {
        '@desc': 'van der Waals radius',
        '@examples': ['vanderwaals >2'],
        regex: reFloat, map: x => parseFloat(x),
        level: 'atom-test', property: B.acp('vdw')
    },
    vectorScale: {
        '@desc': 'vibration vector scale',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    volume: {
        '@desc': 'approximate van der Waals volume for this atom. Note, {*}.volume gives an average; use {*}.volume.sum to get total volume.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    vXyz: {
        '@desc': 'vibration vector, or individual components as %vx %vy %vz. For atoms without vibration vectors, {xx}.vXyz is -1; in a label, %[vxyz] is blank.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    vX: {
        '@desc': 'vibration vector X coordinate; for atoms without vibration vector, {xx}.vX is NaN (same for vY and vZ)',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    vY: {
        '@desc': 'vibration vector Y coordinate',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    vZ: {
        '@desc': 'vibration vector Z coordinate',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
    xyz: {
        '@desc': 'Cartesian XYZ coordinates; select xyz > 1.0 selects atoms more than one Angstrom from the origin.',
        '@examples': [''],
        isUnsupported: true,
        regex: /[0-9\s{}-]+/, map: str,
        level: 'atom-test'
    },
};

