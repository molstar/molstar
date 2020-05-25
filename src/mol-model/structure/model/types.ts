/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import BitFlags from '../../../mol-util/bit-flags';
import { SaccharideCompIdMap } from '../structure/carbohydrates/constants';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { SetUtils } from '../../../mol-util/set';
import { EntitySubtype, ChemicalComponent } from './properties/common';
import { LipidNames } from './types/lipids';

const _esCache = (function () {
    const cache = Object.create(null);
    const letters: string[] = [];
    for (let i = 'A'.charCodeAt(0); i <= 'Z'.charCodeAt(0); i++) letters[letters.length] = String.fromCharCode(i);
    for (let i = 'a'.charCodeAt(0); i <= 'z'.charCodeAt(0); i++) letters[letters.length] = String.fromCharCode(i);
    for (let i = '0'.charCodeAt(0); i <= '9'.charCodeAt(0); i++) letters[letters.length] = String.fromCharCode(i);

    for (const k of letters) {
        cache[k] = k.toUpperCase();
        for (const l of letters) {
            cache[k + l] = (k + l).toUpperCase();
            for (const m of letters) {
                cache[k + l + m] = (k + l + m).toUpperCase();
            }
        }
    }
    return cache;
}());
export type ElementSymbol = string & { '@type': 'element-symbol' }
export function ElementSymbol(s: string): ElementSymbol {
    return _esCache[s] || s.toUpperCase();
}

const _elementByAtomicNumber = new Map(
    ([[1, 'H'], [2, 'He'], [3, 'Li'], [4, 'Be'], [5, 'B'], [6, 'C'], [7, 'N'], [8, 'O'], [9, 'F'], [10, 'Ne'], [11, 'Na'], [12, 'Mg'], [13, 'Al'], [14, 'Si'], [15, 'P'], [16, 'S'], [17, 'Cl'], [18, 'Ar'], [19, 'K'], [20, 'Ca'], [21, 'Sc'], [22, 'Ti'], [23, 'V'], [24, 'Cr'], [25, 'Mn'], [26, 'Fe'], [27, 'Co'], [28, 'Ni'], [29, 'Cu'], [30, 'Zn'], [31, 'Ga'], [32, 'Ge'], [33, 'As'], [34, 'Se'], [35, 'Br'], [36, 'Kr'], [37, 'Rb'], [38, 'Sr'], [39, 'Y'], [40, 'Zr'], [41, 'Nb'], [42, 'Mo'], [43, 'Tc'], [44, 'Ru'], [45, 'Rh'], [46, 'Pd'], [47, 'Ag'], [48, 'Cd'], [49, 'In'], [50, 'Sn'], [51, 'Sb'], [52, 'Te'], [53, 'I'], [54, 'Xe'], [55, 'Cs'], [56, 'Ba'], [57, 'La'], [58, 'Ce'], [59, 'Pr'], [60, 'Nd'], [61, 'Pm'], [62, 'Sm'], [63, 'Eu'], [64, 'Gd'], [65, 'Tb'], [66, 'Dy'], [67, 'Ho'], [68, 'Er'], [69, 'Tm'], [70, 'Yb'], [71, 'Lu'], [72, 'Hf'], [73, 'Ta'], [74, 'W'], [75, 'Re'], [76, 'Os'], [77, 'Ir'], [78, 'Pt'], [79, 'Au'], [80, 'Hg'], [81, 'Tl'], [82, 'Pb'], [83, 'Bi'], [84, 'Po'], [85, 'At'], [86, 'Rn'], [87, 'Fr'], [88, 'Ra'], [89, 'Ac'], [90, 'Th'], [91, 'Pa'], [92, 'U'], [93, 'Np'], [94, 'Pu'], [95, 'Am'], [96, 'Cm'], [97, 'Bk'], [98, 'Cf'], [99, 'Es'], [100, 'Fm'], [101, 'Md'], [102, 'No'], [103, 'Lr'], [104, 'Rf'], [105, 'Db'], [106, 'Sg'], [107, 'Bh'], [108, 'Hs'], [109, 'Mt'], [110, 'Ds'], [111, 'Rg'], [112, 'Cn'], [113, 'Uut'], [114, 'Fl'], [115, 'Uup'], [116, 'Lv'], [117, 'Uus'], [118, 'Uuo']] as const)
        .map(e => [e[0], ElementSymbol(e[1])]));

export function getElementFromAtomicNumber(n: number) {
    if (_elementByAtomicNumber.has(n as any)) return _elementByAtomicNumber.get(n as any)!;
    return ElementSymbol('H');
}

/** Entity types as defined in the mmCIF dictionary */
export const enum EntityType {
    'unknown', 'polymer', 'non-polymer', 'macrolide', 'water', 'branched'
}

export const enum MoleculeType {
    /** The molecule type is not known */
    Unknown,
    /** A known, but here not listed molecule type */
    Other,
    /** Water molecule */
    Water,
    /** Small ionic molecule */
    Ion,
    /** Lipid molecule */
    Lipid,
    /** Protein, e.g. component type included in `ProteinComponentTypeNames` */
    Protein,
    /** RNA, e.g. component type included in `RNAComponentTypeNames` */
    RNA,
    /** DNA, e.g. component type included in `DNAComponentTypeNames` */
    DNA,
    /** PNA, peptide nucleic acid, comp id included in `PeptideBaseNames` */
    PNA,
    /** Saccharide, e.g. component type included in `SaccharideComponentTypeNames` */
    Saccharide
}

export const enum PolymerType {
    /** not applicable */
    NA,
    Protein,
    GammaProtein,
    BetaProtein,
    RNA,
    DNA,
    PNA,
}

export type AtomRole = 'trace' | 'directionFrom' | 'directionTo' | 'backboneStart' | 'backboneEnd' | 'coarseBackbone'

export const PolymerTypeAtomRoleId: { [k in PolymerType]: { [k in AtomRole]: Set<string> } } = {
    [PolymerType.NA]: {
        trace: new Set(),
        directionFrom: new Set(),
        directionTo: new Set(),
        backboneStart: new Set(),
        backboneEnd: new Set(),
        coarseBackbone: new Set()
    },
    [PolymerType.Protein]: {
        trace: new Set(['CA']),
        directionFrom: new Set(['C']),
        directionTo: new Set(['O', 'OC1', 'O1', 'OX1', 'OXT']),
        backboneStart: new Set(['N']),
        backboneEnd: new Set(['C']),
        // CA1 is used e.g. in GFP chromophores
        // BB is often used for coarse grained models
        coarseBackbone: new Set(['CA', 'BB', 'CA1'])
    },
    [PolymerType.GammaProtein]: {
        trace: new Set(['CA']),
        directionFrom: new Set(['C']),
        directionTo: new Set(['O']),
        backboneStart: new Set(['N']),
        backboneEnd: new Set(['CD']),
        coarseBackbone: new Set(['CA'])
    },
    [PolymerType.BetaProtein]: {
        trace: new Set(['CA']),
        directionFrom: new Set(['C']),
        directionTo: new Set(['O']),
        backboneStart: new Set(['N']),
        backboneEnd: new Set(['CG']),
        coarseBackbone: new Set(['CA'])
    },
    [PolymerType.RNA]: {
        trace: new Set(['O3\'', 'O3*']),
        directionFrom: new Set(['C4\'', 'C4*']),
        directionTo: new Set(['C3\'', 'C3*']),
        backboneStart: new Set(['P']),
        backboneEnd: new Set(['O3\'', 'O3*']),
        coarseBackbone: new Set(['P'])
    },
    [PolymerType.DNA]: {
        trace: new Set(['O3\'', 'O3*']),
        directionFrom: new Set(['C3\'', 'C3*']),
        directionTo: new Set(['C1\'', 'C1*']),
        backboneStart: new Set(['P']),
        backboneEnd: new Set(['O3\'', 'O3*']),
        coarseBackbone: new Set(['P'])
    },
    [PolymerType.PNA]: {
        trace: new Set(['N4\'', 'N4*']),
        directionFrom: new Set(['N4\'', 'N4*']),
        directionTo: new Set(['C7\'', 'C7*']),
        backboneStart: new Set(['N1\'', 'N1*']),
        backboneEnd: new Set(['C\'', 'C*']),
        coarseBackbone: new Set(['P'])
    }
};

export const ProteinBackboneAtoms = new Set([
    'CA', 'C', 'N', 'O',
    'O1', 'O2', 'OC1', 'OC2', 'OT1', 'OT2', 'OX1', 'OXT',
    'H', 'H1', 'H2', 'H3', 'HA', 'HN', 'HXT',
    'BB'
]);

export const NucleicBackboneAtoms = new Set([
    'P', 'OP1', 'OP2', 'HOP2', 'HOP3',
    'O2\'', 'O3\'', 'O4\'', 'O5\'', 'C1\'', 'C2\'', 'C3\'', 'C4\'', 'C5\'',
    'H1\'', 'H2\'', 'H2\'\'', 'HO2\'', 'H3\'', 'H4\'', 'H5\'', 'H5\'\'', 'HO3\'', 'HO5\'',
    'O2*', 'O3*', 'O4*', 'O5*', 'C1*', 'C2*', 'C3*', 'C4*', 'C5*'
]);

/** Chemical component type names for D-linked protein */
export const DProteinComponentTypeNames = new Set([
    'D-PEPTIDE LINKING', 'D-PEPTIDE NH3 AMINO TERMINUS',
    'D-PEPTIDE COOH CARBOXY TERMINUS', 'D-GAMMA-PEPTIDE, C-DELTA LINKING',
    'D-BETA-PEPTIDE, C-GAMMA LINKING'
]);

/** Chemical component type names for L-linked protein */
export const LProteinComponentTypeNames = new Set([
    'L-PEPTIDE LINKING', 'L-PEPTIDE NH3 AMINO TERMINUS',
    'L-PEPTIDE COOH CARBOXY TERMINUS', 'L-GAMMA-PEPTIDE, C-DELTA LINKING',
    'L-BETA-PEPTIDE, C-GAMMA LINKING'
]);

/** Chemical component type names for gamma protein, overlaps with D/L-linked */
export const GammaProteinComponentTypeNames = new Set([
    'D-GAMMA-PEPTIDE, C-DELTA LINKING', 'L-GAMMA-PEPTIDE, C-DELTA LINKING'
]);

/** Chemical component type names for beta protein, overlaps with D/L-linked */
export const BetaProteinComponentTypeNames = new Set([
    'D-BETA-PEPTIDE, C-GAMMA LINKING', 'L-BETA-PEPTIDE, C-GAMMA LINKING'
]);

/** Chemical component type names for protein termini, overlaps with D/L-linked */
export const ProteinTerminusComponentTypeNames = new Set([
    'D-PEPTIDE NH3 AMINO TERMINUS', 'D-PEPTIDE COOH CARBOXY TERMINUS',
    'L-PEPTIDE NH3 AMINO TERMINUS', 'L-PEPTIDE COOH CARBOXY TERMINUS'
]);

/** Chemical component type names for pepdite-like protein */
export const OtherProteinComponentTypeNames = new Set([
    'PEPTIDE LINKING', 'PEPTIDE-LIKE',
]);

/** Chemical component type names for protein */
export const ProteinComponentTypeNames = SetUtils.unionMany(
    DProteinComponentTypeNames, LProteinComponentTypeNames, OtherProteinComponentTypeNames
);

/** Chemical component type names for DNA */
export const DNAComponentTypeNames = new Set([
    'DNA LINKING', 'L-DNA LINKING', 'DNA OH 5 PRIME TERMINUS', 'DNA OH 3 PRIME TERMINUS',
]);

/** Chemical component type names for RNA */
export const RNAComponentTypeNames = new Set([
    'RNA LINKING', 'L-RNA LINKING', 'RNA OH 5 PRIME TERMINUS', 'RNA OH 3 PRIME TERMINUS',
]);

/** Chemical component type names for saccharide */
export const SaccharideComponentTypeNames = new Set([
    'D-SACCHARIDE, BETA LINKING', 'L-SACCHARIDE, BETA LINKING',
    'D-SACCHARIDE, ALPHA LINKING', 'L-SACCHARIDE, ALPHA LINKING',
    'L-SACCHARIDE', 'D-SACCHARIDE', 'SACCHARIDE',
    // the following four are marked to be deprecated in the mmCIF dictionary
    'D-SACCHARIDE 1,4 AND 1,4 LINKING', 'L-SACCHARIDE 1,4 AND 1,4 LINKING',
    'D-SACCHARIDE 1,4 AND 1,6 LINKING', 'L-SACCHARIDE 1,4 AND 1,6 LINKING',
]);

/** Chemical component type names for other */
export const OtherComponentTypeNames = new Set([
    'NON-POLYMER', 'OTHER'
]);

/** Common names for water molecules */
export const WaterNames = new Set([
    'SOL', 'WAT', 'HOH', 'H2O', 'W', 'DOD', 'D3O', 'TIP', 'TIP3', 'TIP4', 'SPC'
]);

export const AminoAcidNamesL = new Set([
    'HIS', 'ARG', 'LYS', 'ILE', 'PHE', 'LEU', 'TRP', 'ALA', 'MET', 'PRO', 'CYS',
    'ASN', 'VAL', 'GLY', 'SER', 'GLN', 'TYR', 'ASP', 'GLU', 'THR', 'SEC', 'PYL',
    'UNK' // unknown amino acid from CCD
]);
export const AminoAcidNamesD = new Set([
    'DAL', // D-ALANINE
    'DAR', // D-ARGININE
    'DSG', // D-ASPARAGINE
    'DAS', // D-ASPARTIC ACID
    'DCY', // D-CYSTEINE
    'DGL', // D-GLUTAMIC ACID
    'DGN', // D-GLUTAMINE
    'DHI', // D-HISTIDINE
    'DIL', // D-ISOLEUCINE
    'DLE', // D-LEUCINE
    'DLY', // D-LYSINE
    'MED', // D-METHIONINE
    'DPN', // D-PHENYLALANINE
    'DPR', // D-PROLINE
    'DSN', // D-SERINE
    'DTH', // D-THREONINE
    'DTR', // D-TRYPTOPHAN
    'DTY', // D-TYROSINE
    'DVA', // D-VALINE
    'DNE' // D-NORLEUCINE
    // ???  // D-SELENOCYSTEINE
]);
export const AminoAcidNames = SetUtils.unionMany(AminoAcidNamesL, AminoAcidNamesD);

export const RnaBaseNames = new Set([
    'A', 'C', 'T', 'G', 'I', 'U',
    'N' // unknown RNA base from CCD
]);
export const DnaBaseNames = new Set([
    'DA', 'DC', 'DT', 'DG', 'DI', 'DU',
    'DN' // unknown DNA base from CCD
]);
export const PeptideBaseNames = new Set(['APN', 'CPN', 'TPN', 'GPN']);
export const PurineBaseNames = new Set(['A', 'G', 'I', 'DA', 'DG', 'DI', 'APN', 'GPN']);
export const PyrimidineBaseNames = new Set(['C', 'T', 'U', 'DC', 'DT', 'DU', 'CPN', 'TPN']);
export const BaseNames = SetUtils.unionMany(RnaBaseNames, DnaBaseNames, PeptideBaseNames);

export const isPurineBase = (compId: string) => PurineBaseNames.has(compId.toUpperCase());
export const isPyrimidineBase = (compId: string) => PyrimidineBaseNames.has(compId.toUpperCase());

export const PolymerNames = SetUtils.unionMany(AminoAcidNames, BaseNames);

/** get the molecule type from component type and id */
export function getMoleculeType(compType: string, compId: string): MoleculeType {
    compType = compType.toUpperCase();
    compId = compId.toUpperCase();
    if (PeptideBaseNames.has(compId)) {
        return MoleculeType.PNA;
    } else if (ProteinComponentTypeNames.has(compType)) {
        return MoleculeType.Protein;
    } else if (RNAComponentTypeNames.has(compType)) {
        return MoleculeType.RNA;
    } else if (DNAComponentTypeNames.has(compType)) {
        return MoleculeType.DNA;
    } else if (SaccharideComponentTypeNames.has(compType)) {
        return MoleculeType.Saccharide;
    } else if (WaterNames.has(compId)) {
        return MoleculeType.Water;
    } else if (IonNames.has(compId)) {
        return MoleculeType.Ion;
    } else if (LipidNames.has(compId)) {
        return MoleculeType.Lipid;
    } else if (OtherComponentTypeNames.has(compType)) {
        if (SaccharideCompIdMap.has(compId)) {
            // trust our saccharide table more than given 'non-polymer' or 'other' component type
            return MoleculeType.Saccharide;
        } else {
            return MoleculeType.Other;
        }
    } else {
        return MoleculeType.Unknown;
    }
}

export function getPolymerType(compType: string, molType: MoleculeType): PolymerType {
    compType = compType.toUpperCase();
    if (molType === MoleculeType.Protein) {
        if (GammaProteinComponentTypeNames.has(compType)) {
            return PolymerType.GammaProtein;
        } else if (BetaProteinComponentTypeNames.has(compType)) {
            return PolymerType.BetaProtein;
        } else if (ProteinTerminusComponentTypeNames.has(compType)) {
            return PolymerType.NA;
        } else {
            return PolymerType.Protein;
        }
    } else if (molType === MoleculeType.RNA) {
        return PolymerType.RNA;
    } else if (molType === MoleculeType.DNA) {
        return PolymerType.DNA;
    } else if (molType === MoleculeType.PNA) {
        return PolymerType.PNA;
    } else {
        return PolymerType.NA;
    }
}

export function getComponentType(compId: string): mmCIF_Schema['chem_comp']['type']['T'] {
    compId = compId.toUpperCase();
    if (AminoAcidNames.has(compId)) {
        return 'peptide linking';
    } else if (RnaBaseNames.has(compId)) {
        return 'RNA linking';
    } else if (DnaBaseNames.has(compId)) {
        return 'DNA linking';
    } else if (SaccharideCompIdMap.has(compId)) {
        return 'saccharide';
    } else {
        return 'other';
    }
}

export function getDefaultChemicalComponent(compId: string): ChemicalComponent {
    // TODO: this is to make the chem_comp_type property work if chem_comp category is not present.
    // should we try to set the formula etc better?
    return {
        formula: '',
        formula_weight: 0,
        id: compId,
        name: compId,
        mon_nstd_flag: PolymerNames.has(compId) ? 'y' : 'n',
        pdbx_synonyms: [],
        type: getComponentType(compId)
    };
}

export function getEntityType(compId: string): mmCIF_Schema['entity']['type']['T'] {
    compId = compId.toUpperCase();
    if (WaterNames.has(compId)) {
        return 'water';
    } else if (PolymerNames.has(compId)) {
        return 'polymer';
    } else if (SaccharideCompIdMap.has(compId)) {
        return 'branched';
    } else {
        return 'non-polymer';
    }
}

export function getEntitySubtype(compId: string, compType: string): EntitySubtype {
    compId = compId.toUpperCase();
    compType = compType.toUpperCase();
    if (LProteinComponentTypeNames.has(compType)) {
        return 'polypeptide(L)';
    } else if (DProteinComponentTypeNames.has(compType)) {
        return 'polypeptide(D)';
    } else if (RNAComponentTypeNames.has(compType)) {
        return 'polyribonucleotide';
    } else if (DNAComponentTypeNames.has(compType)) {
        return 'polydeoxyribonucleotide';
    } else if (SaccharideComponentTypeNames.has(compType)) {
        return 'oligosaccharide';
    } else if (SaccharideCompIdMap.has(compId)) {
        return 'oligosaccharide';
    } else if (PeptideBaseNames.has(compId)) {
        return 'peptide nucleic acid';
    } else if (AminoAcidNamesL.has(compId)) {
        return 'polypeptide(L)';
    } else if (AminoAcidNamesD.has(compId)) {
        return 'polypeptide(D)';
    } else if (RnaBaseNames.has(compId)) {
        return 'polyribonucleotide';
    } else if (DnaBaseNames.has(compId)) {
        return 'polydeoxyribonucleotide';
    } else if (IonNames.has(compId)) {
        return 'ion';
    } else if (LipidNames.has(compId)) {
        return 'lipid';
    } else {
        return 'other';
    }
}

export function isPolymer(moleculeType: MoleculeType) {
    return isNucleic(moleculeType) || isProtein(moleculeType);
}

export function isNucleic(moleculeType: MoleculeType) {
    return moleculeType === MoleculeType.DNA || moleculeType === MoleculeType.RNA || moleculeType === MoleculeType.PNA;
}

export function isProtein(moleculeType: MoleculeType) {
    return moleculeType === MoleculeType.Protein;
}

/**
 * TODO write script that read CCD and outputs list of ion names
 *
 * all chemical components with the word "ion" in their name, Sep 2016
 *
 * SET SESSION group_concat_max_len = 1000000;
 * SELECT GROUP_CONCAT(id_ ORDER BY id_ ASC SEPARATOR '", "') from
 * (
 *     SELECT count(obj_id) as c, id_
 *     FROM pdb.chem_comp WHERE name LIKE "% ION%"
 *     GROUP BY id_
 * ) AS t1;
 */
export const IonNames = new Set([
    '118', '119', '1AL', '1CU', '2FK', '2HP', '2OF', '3CO',
    '3MT', '3NI', '3OF', '3P8', '4MO', '4PU', '543', '6MO', 'ACT', 'AG', 'AL',
    'ALF', 'AM', 'ATH', 'AU', 'AU3', 'AUC', 'AZI', 'BA', 'BCT', 'BEF', 'BF4', 'BO4',
    'BR', 'BS3', 'BSY', 'CA', 'CAC', 'CD', 'CD1', 'CD3', 'CD5', 'CE', 'CHT', 'CL',
    'CO', 'CO3', 'CO5', 'CON', 'CR', 'CS', 'CSB', 'CU', 'CU1', 'CU3', 'CUA', 'CUZ',
    'CYN', 'DME', 'DMI', 'DSC', 'DTI', 'DY', 'E4N', 'EDR', 'EMC', 'ER3', 'EU',
    'EU3', 'F', 'FE', 'FE2', 'FPO', 'GA', 'GD3', 'GEP', 'HAI', 'HG', 'HGC', 'IN',
    'IOD', 'IR', 'IR3', 'IRI', 'IUM', 'K', 'KO4', 'LA', 'LCO', 'LCP', 'LI', 'LU',
    'MAC', 'MG', 'MH2', 'MH3', 'MLI', 'MLT', 'MMC', 'MN', 'MN3', 'MN5', 'MN6',
    'MO1', 'MO2', 'MO3', 'MO4', 'MO5', 'MO6', 'MOO', 'MOS', 'MOW', 'MW1', 'MW2',
    'MW3', 'NA', 'NA2', 'NA5', 'NA6', 'NAO', 'NAW', 'NCO', 'NET', 'NH4', 'NI',
    'NI1', 'NI2', 'NI3', 'NO2', 'NO3', 'NRU', 'O4M', 'OAA', 'OC1', 'OC2', 'OC3',
    'OC4', 'OC5', 'OC6', 'OC7', 'OC8', 'OCL', 'OCM', 'OCN', 'OCO', 'OF1', 'OF2',
    'OF3', 'OH', 'OS', 'OS4', 'OXL', 'PB', 'PBM', 'PD', 'PDV', 'PER', 'PI', 'PO3',
    'PO4', 'PR', 'PT', 'PT4', 'PTN', 'RB', 'RH3', 'RHD', 'RU', 'SB', 'SCN', 'SE4',
    'SEK', 'SM', 'SMO', 'SO3', 'SO4', 'SR', 'T1A', 'TB', 'TBA', 'TCN', 'TEA', 'TH',
    'THE', 'TL', 'TMA', 'TRA', 'UNX', 'V', 'VN3', 'VO4', 'W', 'WO5', 'Y1', 'YB',
    'YB2', 'YH', 'YT3', 'ZCM', 'ZN', 'ZN2', 'ZN3', 'ZNO', 'ZO3',
    // additional ion names
    'OHX'
]);

export type SecondaryStructureType = BitFlags<SecondaryStructureType.Flag>
export namespace SecondaryStructureType {
    export const is: (ss: SecondaryStructureType, f: Flag) => boolean = BitFlags.has;
    export const create: (fs: Flag) => SecondaryStructureType = BitFlags.create;

    export const enum Flag {
        None = 0x0,

        // category
        DoubleHelix = 0x1,
        Helix = 0x2,
        Beta = 0x4,
        Bend = 0x8,
        Turn = 0x10,

        // category variant
        LeftHanded = 0x20,  // helix
        RightHanded = 0x40,

        ClassicTurn = 0x80,  // turn
        InverseTurn = 0x100,

        // sub-category
        HelixOther = 0x200,  // protein
        Helix27 = 0x400,
        Helix3Ten = 0x800,
        HelixAlpha = 0x1000,
        HelixGamma = 0x2000,
        HelixOmega = 0x4000,
        HelixPi = 0x8000,
        HelixPolyproline = 0x10000,

        DoubleHelixOther = 0x20000,  // nucleic
        DoubleHelixZ = 0x40000,
        DoubleHelixA = 0x80000,
        DoubleHelixB = 0x100000,

        BetaOther = 0x200000,  // protein
        BetaStrand = 0x400000,  // single strand
        BetaSheet = 0x800000,  // multiple hydrogen bonded strands
        BetaBarell = 0x1000000,  // closed series of sheets

        TurnOther = 0x2000000,  // protein
        Turn1 = 0x4000000,
        Turn2 = 0x8000000,
        Turn3 = 0x10000000,

        NA = 0x20000000,  // not applicable/available
    }

    export const SecondaryStructureMmcif: { [value: string]: number } = {
        HELX_LH_27_P: Flag.Helix | Flag.LeftHanded | Flag.Helix27,  // left-handed 2-7 helix (protein)
        HELX_LH_3T_P: Flag.Helix | Flag.LeftHanded | Flag.Helix3Ten,  // left-handed 3-10 helix (protein)
        HELX_LH_AL_P: Flag.Helix | Flag.LeftHanded | Flag.HelixAlpha,  // left-handed alpha helix (protein)
        HELX_LH_A_N: Flag.DoubleHelix | Flag.LeftHanded | Flag.DoubleHelixA,  // left-handed A helix (nucleic acid)
        HELX_LH_B_N: Flag.DoubleHelix | Flag.LeftHanded | Flag.DoubleHelixB,  // left-handed B helix (nucleic acid)
        HELX_LH_GA_P: Flag.Helix | Flag.LeftHanded | Flag.HelixGamma,  // left-handed gamma helix (protein)
        HELX_LH_N: Flag.DoubleHelix | Flag.LeftHanded,  // left-handed helix with type not specified (nucleic acid)
        HELX_LH_OM_P: Flag.Helix | Flag.LeftHanded | Flag.HelixOmega,  // left-handed omega helix (protein)
        HELX_LH_OT_N: Flag.DoubleHelix | Flag.LeftHanded | Flag.DoubleHelixOther,  // left-handed helix with type that does not conform to an accepted category (nucleic acid)
        HELX_LH_OT_P: Flag.Helix | Flag.LeftHanded | Flag.HelixOther,  // left-handed helix with type that does not conform to an accepted category (protein)
        HELX_LH_P: Flag.Helix | Flag.LeftHanded,  // left-handed helix with type not specified (protein)
        HELX_LH_PI_P: Flag.Helix | Flag.LeftHanded | Flag.HelixPi,  // left-handed pi helix (protein)
        HELX_LH_PP_P: Flag.Helix | Flag.LeftHanded | Flag.HelixPolyproline,  // left-handed polyproline helix (protein)
        HELX_LH_Z_N: Flag.DoubleHelix | Flag.LeftHanded | Flag.DoubleHelixZ,  // left-handed Z helix (nucleic acid)
        HELX_N: Flag.DoubleHelix,  // helix with handedness and type not specified (nucleic acid)
        HELX_OT_N: Flag.DoubleHelix,  // helix with handedness and type that do not conform to an accepted category (nucleic acid)
        HELX_OT_P: Flag.Helix,  // helix with handedness and type that do not conform to an accepted category (protein)
        HELX_P: Flag.Helix,  // helix with handedness and type not specified (protein)
        HELX_RH_27_P: Flag.Helix | Flag.RightHanded | Flag.Helix27,  // right-handed 2-7 helix (protein)
        HELX_RH_3T_P: Flag.Helix | Flag.RightHanded | Flag.Helix3Ten,  // right-handed 3-10 helix (protein)
        HELX_RH_AL_P: Flag.Helix | Flag.RightHanded | Flag.HelixAlpha,  // right-handed alpha helix (protein)
        HELX_RH_A_N: Flag.DoubleHelix | Flag.RightHanded | Flag.DoubleHelixA,  // right-handed A helix (nucleic acid)
        HELX_RH_B_N: Flag.DoubleHelix | Flag.RightHanded | Flag.DoubleHelixB,  // right-handed B helix (nucleic acid)
        HELX_RH_GA_P: Flag.Helix | Flag.RightHanded | Flag.HelixGamma,  // right-handed gamma helix (protein)
        HELX_RH_N: Flag.DoubleHelix | Flag.RightHanded,  // right-handed helix with type not specified (nucleic acid)
        HELX_RH_OM_P: Flag.Helix | Flag.RightHanded | Flag.HelixOmega,  // right-handed omega helix (protein)
        HELX_RH_OT_N: Flag.DoubleHelix | Flag.RightHanded | Flag.DoubleHelixOther,  // right-handed helix with type that does not conform to an accepted category (nucleic acid)
        HELX_RH_OT_P: Flag.Helix | Flag.RightHanded | Flag.HelixOther,  // right-handed helix with type that does not conform to an accepted category (protein)
        HELX_RH_P: Flag.Helix | Flag.RightHanded,  // right-handed helix with type not specified (protein)
        HELX_RH_PI_P: Flag.Helix | Flag.RightHanded | Flag.HelixPi,  // right-handed pi helix (protein)
        HELX_RH_PP_P: Flag.Helix | Flag.RightHanded | Flag.HelixPolyproline,  // right-handed polyproline helix (protein)
        HELX_RH_Z_N: Flag.DoubleHelix | Flag.RightHanded | Flag.DoubleHelixZ,  // right-handed Z helix (nucleic acid)
        STRN: Flag.Beta | Flag.BetaStrand,  // beta strand (protein)
        TURN_OT_P: Flag.Turn | Flag.TurnOther,  // turn with type that does not conform to an accepted category (protein)
        TURN_P: Flag.Turn,  // turn with type not specified (protein)
        TURN_TY1P_P: Flag.Turn | Flag.InverseTurn | Flag.Turn1,  // type I prime turn (protein)
        TURN_TY1_P: Flag.Turn | Flag.ClassicTurn | Flag.Turn1,  // type I turn (protein)
        TURN_TY2P_P: Flag.Turn | Flag.InverseTurn | Flag.Turn2,  // type II prime turn (protein)
        TURN_TY2_P: Flag.Turn | Flag.ClassicTurn | Flag.Turn2,  // type II turn (protein)
        TURN_TY3P_P: Flag.Turn | Flag.InverseTurn | Flag.Turn3,  // type III prime turn (protein)
        TURN_TY3_P: Flag.Turn | Flag.ClassicTurn | Flag.Turn3,  // type III turn (protein)
    };

    export const SecondaryStructurePdb: { [value: string]: number } = {
        1: Flag.Helix | Flag.RightHanded | Flag.HelixAlpha,  // Right-handed alpha (default)
        2: Flag.Helix | Flag.RightHanded | Flag.HelixOmega,  // Right-handed omega
        3: Flag.Helix | Flag.RightHanded | Flag.HelixPi,  // Right-handed pi
        4: Flag.Helix | Flag.RightHanded | Flag.HelixGamma,  // Right-handed gamma
        5: Flag.Helix | Flag.RightHanded | Flag.Helix3Ten,  // Right-handed 310
        6: Flag.Helix | Flag.LeftHanded | Flag.HelixAlpha,  // Left-handed alpha
        7: Flag.Helix | Flag.LeftHanded | Flag.HelixOmega,  // Left-handed omega
        8: Flag.Helix | Flag.LeftHanded | Flag.HelixGamma,  // Left-handed gamma
        9: Flag.Helix | Flag.Helix27,  // 27 ribbon/helix
        10: Flag.Helix | Flag.HelixPolyproline,  // Polyproline
    };

    export const SecondaryStructureStride: { [value: string]: number } = {
        H: Flag.Helix | Flag.HelixAlpha,  // Alpha helix
        G: Flag.Helix | Flag.Helix3Ten,  // 3-10 helix
        I: Flag.Helix | Flag.HelixPi,  // PI-helix
        E: Flag.Beta | Flag.BetaSheet,  // Extended conformation
        B: Flag.Beta | Flag.BetaStrand,  // Isolated bridge
        T: Flag.Turn,  // Turn
        C: Flag.NA,  // Coil (none of the above)
    };

    export const SecondaryStructureDssp: { [value: string]: number } = {
        H: Flag.Helix | Flag.HelixAlpha,  // alpha-helix
        B: Flag.Beta | Flag.BetaStrand,  // residue in isolated beta-bridge
        E: Flag.Beta | Flag.BetaSheet,  // extended strand, participates in beta ladder
        G: Flag.Helix | Flag.Helix3Ten,  // 3-helix (310 helix)
        I: Flag.Helix | Flag.HelixPi,  // 5 helix (pi-helix)
        T: Flag.Turn,  // hydrogen bonded turn
        S: Flag.Bend,  // bend
    };
}

/** Maximum accessible surface area observed for amino acids. Taken from: http://dx.doi.org/10.1371/journal.pone.0080635 */
export const MaxAsa = {
    'ALA': 121.0,
    'ARG': 265.0,
    'ASN': 187.0,
    'ASP': 187.0,
    'CYS': 148.0,
    'GLU': 214.0,
    'GLN': 214.0,
    'GLY': 97.0,
    'HIS': 216.0,
    'ILE': 195.0,
    'LEU': 191.0,
    'LYS': 230.0,
    'MET': 203.0,
    'PHE': 228.0,
    'PRO': 154.0,
    'SER': 143.0,
    'THR': 163.0,
    'TRP': 264.0,
    'TYR': 255.0,
    'VAL': 165.0
};
export const DefaultMaxAsa = 121.0;

export type BondType = BitFlags<BondType.Flag>
export namespace BondType {
    export const is: (b: BondType, f: Flag) => boolean = BitFlags.has;
    export const enum Flag {
        None = 0x0,
        Covalent = 0x1,
        MetallicCoordination = 0x2,
        HydrogenBond = 0x4,
        Disulfide = 0x8,
        Aromatic = 0x10,
        Computed = 0x20
        // currently at most 16 flags are supported!!
    }

    export function create(flags: Flag): BondType {
        return BitFlags.create(flags);
    }

    export function isCovalent(flags: BondType.Flag) {
        return (flags & BondType.Flag.Covalent) !== 0;
    }

    export const Names = {
        'covalent': Flag.Covalent,
        'metal-coordination': Flag.MetallicCoordination,
        'hydrogen-bond': Flag.HydrogenBond,
        'disulfide': Flag.HydrogenBond,
        'aromatic': Flag.HydrogenBond,
        'computed': Flag.HydrogenBond,
    };
    export type Names = keyof typeof Names

    export function isName(name: string): name is Names {
        return name in Names;
    }

    export function fromName(name: Names): Flag {
        switch (name) {
            case 'covalent': return Flag.Covalent;
            case 'metal-coordination': return Flag.MetallicCoordination;
            case 'hydrogen-bond': return Flag.HydrogenBond;
            case 'disulfide': return Flag.Disulfide;
            case 'aromatic': return Flag.Aromatic;
            case 'computed': return Flag.Computed;
        }
    }

    export function fromNames(names: Names[]): Flag {
        let f = Flag.None;
        for (let i = 0, il = names.length; i < il; ++i) {
            f |= fromName(names[i]);
        }
        return f;
    }
}

/**
 * "Experimentally determined hydrophobicity scale for proteins at membrane interfaces"
 * by Wimely and White (doi:10.1038/nsb1096-842)
 * http://blanco.biomol.uci.edu/Whole_residue_HFscales.txt
 * https://www.nature.com/articles/nsb1096-842
 */
export const ResidueHydrophobicity = {
    // AA  DGwif   DGwoct  Oct-IF
    'ALA': [0.17, 0.50, 0.33],
    'ARG': [0.81, 1.81, 1.00],
    'ASN': [0.42, 0.85, 0.43],
    'ASP': [1.23, 3.64, 2.41],
    'ASH': [-0.07, 0.43, 0.50],
    'CYS': [-0.24, -0.02, 0.22],
    'GLN': [0.58, 0.77, 0.19],
    'GLU': [2.02, 3.63, 1.61],
    'GLH': [-0.01, 0.11, 0.12],
    'GLY': [0.01, 1.15, 1.14],
    // "His+": [  0.96,  2.33,  1.37 ],
    'HIS': [0.17, 0.11, -0.06],
    'ILE': [-0.31, -1.12, -0.81],
    'LEU': [-0.56, -1.25, -0.69],
    'LYS': [0.99, 2.80, 1.81],
    'MET': [-0.23, -0.67, -0.44],
    'PHE': [-1.13, -1.71, -0.58],
    'PRO': [0.45, 0.14, -0.31],
    'SER': [0.13, 0.46, 0.33],
    'THR': [0.14, 0.25, 0.11],
    'TRP': [-1.85, -2.09, -0.24],
    'TYR': [-0.94, -0.71, 0.23],
    'VAL': [0.07, -0.46, -0.53]
};
export const DefaultResidueHydrophobicity = [0.00, 0.00, 0.00];