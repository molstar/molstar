/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import BitFlags from 'mol-util/bit-flags'
import { SaccharideCompIdMap } from '../structure/carbohydrates/constants';
import { mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif';
import { SetUtils } from 'mol-util/set';

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

/** Entity types as defined in the mmCIF dictionary */
export const enum EntityType {
    'unknown', 'polymer', 'non-polymer', 'macrolide', 'water', 'branched'
}

export const enum MoleculeType {
    /** the molecule type is not known */
    unknown,
    /** a known, but here not listed molecule type */
    other,
    /** water molecule */
    water,
    /** small ionic molecule */
    ion,
    /** protein, e.g. component type included in `ProteinComponentTypeNames` */
    protein,
    /** RNA, e.g. component type included in `RNAComponentTypeNames` */
    RNA,
    /** DNA, e.g. component type included in `DNAComponentTypeNames` */
    DNA,
    /** PNA, peptide nucleic acid, comp id included in `PeptideBaseNames` */
    PNA,
    /** sacharide, e.g. component type included in `SaccharideComponentTypeNames` */
    saccharide
}

export type AtomRole = 'trace' | 'direction' | 'backboneStart' | 'backboneEnd' | 'coarseBackbone'

export const MoleculeTypeAtomRoleId: { [k: number]: { [k in AtomRole]: Set<string> } } = {
    [MoleculeType.protein]: {
        trace: new Set(['CA']),
        direction: new Set(['O', 'OC1', 'O1', 'OX1', 'OXT']),
        backboneStart: new Set(['N']),
        backboneEnd: new Set(['C']),
        coarseBackbone: new Set(['CA', 'BB'])
    },
    [MoleculeType.RNA]: {
        trace: new Set(['C4\'', 'C4*']),
        direction: new Set(['C3\'', 'C3*']),
        backboneStart: new Set(['P']),
        backboneEnd: new Set(['O3\'', 'O3*']),
        coarseBackbone: new Set(['P'])
    },
    [MoleculeType.DNA]: {
        trace: new Set(['C3\'', 'C3*']),
        direction: new Set(['C1\'', 'C1*']),
        backboneStart: new Set(['P']),
        backboneEnd: new Set(['O3\'', 'O3*']),
        coarseBackbone: new Set(['P'])
    },
    [MoleculeType.PNA]: {
        trace: new Set(['N4\'', 'N4*']),
        direction: new Set(['C7\'', 'C7*']),
        backboneStart: new Set(['N1\'', 'N1*']),
        backboneEnd: new Set(['C\'', 'C*']),
        coarseBackbone: new Set(['P'])
    }
}

export const ProteinBackboneAtoms = new Set([
    'CA', 'C', 'N', 'O',
    'O1', 'O2', 'OC1', 'OC2', 'OX1', 'OXT',
    'H', 'H1', 'H2', 'H3', 'HA', 'HN',
    'BB'
])

export const NucleicBackboneAtoms = new Set([
    'P', 'OP1', 'OP2', 'HOP2', 'HOP3',
    'O2\'', 'O3\'', 'O4\'', 'O5\'', 'C1\'', 'C2\'', 'C3\'', 'C4\'', 'C5\'',
    'H1\'', 'H2\'', 'H2\'\'', 'HO2\'', 'H3\'', 'H4\'', 'H5\'', 'H5\'\'', 'HO3\'', 'HO5\'',
    'O2*', 'O3*', 'O4*', 'O5*', 'C1*', 'C2*', 'C3*', 'C4*', 'C5*'
])

/** Chemical component type names for protein */
export const ProteinComponentTypeNames = new Set([
    'D-PEPTIDE LINKING', 'L-PEPTIDE LINKING', 'D-PEPTIDE NH3 AMINO TERMINUS',
    'L-PEPTIDE NH3 AMINO TERMINUS', 'D-PEPTIDE COOH CARBOXY TERMINUS',
    'L-PEPTIDE COOH CARBOXY TERMINUS', 'PEPTIDE LINKING', 'PEPTIDE-LIKE',
    'L-GAMMA-PEPTIDE, C-DELTA LINKING', 'D-GAMMA-PEPTIDE, C-DELTA LINKING',
    'L-BETA-PEPTIDE, C-GAMMA LINKING', 'D-BETA-PEPTIDE, C-GAMMA LINKING',
])

/** Chemical component type names for DNA */
export const DNAComponentTypeNames = [
    'DNA LINKING', 'L-DNA LINKING', 'DNA OH 5 PRIME TERMINUS', 'DNA OH 3 PRIME TERMINUS',
]

/** Chemical component type names for RNA */
export const RNAComponentTypeNames = new Set([
    'RNA LINKING', 'L-RNA LINKING', 'RNA OH 5 PRIME TERMINUS', 'RNA OH 3 PRIME TERMINUS',
])

/** Chemical component type names for saccharide */
export const SaccharideComponentTypeNames = new Set([
    'D-SACCHARIDE 1,4 AND 1,4 LINKING', 'L-SACCHARIDE 1,4 AND 1,4 LINKING',
    'D-SACCHARIDE 1,4 AND 1,6 LINKING', 'L-SACCHARIDE 1,4 AND 1,6 LINKING', 'L-SACCHARIDE',
    'D-SACCHARIDE', 'SACCHARIDE',
])

/** Chemical component type names for other */
export const OtherComponentTypeNames = new Set([
    'NON-POLYMER', 'OTHER'
])

/** Common names for water molecules */
export const WaterNames = new Set([
    'SOL', 'WAT', 'HOH', 'H2O', 'W', 'DOD', 'D3O', 'TIP3', 'TIP4', 'SPC'
])

export const AminoAcidNames = new Set([
    'HIS', 'ARG', 'LYS', 'ILE', 'PHE', 'LEU', 'TRP', 'ALA', 'MET', 'PRO', 'CYS',
    'ASN', 'VAL', 'GLY', 'SER', 'GLN', 'TYR', 'ASP', 'GLU', 'THR', 'SEC', 'PYL',

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
])

export const RnaBaseNames = new Set([ 'A', 'C', 'T', 'G', 'I', 'U' ])
export const DnaBaseNames = new Set([ 'DA', 'DC', 'DT', 'DG', 'DI', 'DU' ])
export const PeptideBaseNames = new Set([ 'APN', 'CPN', 'TPN', 'GPN' ])
export const PurinBaseNames = new Set([ 'A', 'G', 'DA', 'DG', 'DI', 'APN', 'GPN' ])
export const PyrimidineBaseNames = new Set([ 'C', 'T', 'U', 'DC', 'DT', 'DU', 'CPN', 'TPN' ])
export const BaseNames = SetUtils.unionMany(RnaBaseNames, DnaBaseNames, PeptideBaseNames)

export const isPurinBase = (compId: string) => PurinBaseNames.has(compId.toUpperCase())
export const isPyrimidineBase = (compId: string) => PyrimidineBaseNames.has(compId.toUpperCase())

/** get the molecule type from component type and id */
export function getMoleculeType(compType: string, compId: string) {
    compType = compType.toUpperCase()
    compId = compId.toUpperCase()
    if (PeptideBaseNames.has(compId)) {
        return MoleculeType.PNA
    } else if (ProteinComponentTypeNames.has(compType)) {
        return MoleculeType.protein
    } else if (RNAComponentTypeNames.has(compType)) {
        return MoleculeType.RNA
    } else if (DNAComponentTypeNames.includes(compType)) {
        return MoleculeType.DNA
    } else if (SaccharideComponentTypeNames.has(compType)) {
        return MoleculeType.saccharide
    } else if (WaterNames.has(compId)) {
        return MoleculeType.water
    } else if (IonNames.has(compId)) {
        return MoleculeType.ion
    } else if (OtherComponentTypeNames.has(compType)) {
        return MoleculeType.other
    } else {
        return MoleculeType.unknown
    }
}

export function getComponentType(compId: string): mmCIF_Schema['chem_comp']['type']['T'] {
    compId = compId.toUpperCase()
    if (AminoAcidNames.has(compId)) {
        return 'peptide linking'
    } else if (RnaBaseNames.has(compId)) {
        return 'RNA linking'
    } else if (DnaBaseNames.has(compId)) {
        return 'DNA linking'
    } else if (SaccharideCompIdMap.has(compId)) {
        return 'saccharide'
    } else {
        return 'other'
    }
}

export function getEntityType(compId: string): mmCIF_Schema['entity']['type']['T'] {
    compId = compId.toUpperCase()
    if (AminoAcidNames.has(compId) || RnaBaseNames.has(compId) || DnaBaseNames.has(compId)) {
        return 'polymer'
    } else if (SaccharideCompIdMap.has(compId)) {
        return 'branched'
    } else if (WaterNames.has(compId)) {
        return 'water'
    } else {
        return 'non-polymer'
    }
}

export function isPolymer(moleculeType: MoleculeType) {
    return moleculeType === MoleculeType.protein || moleculeType === MoleculeType.DNA || moleculeType === MoleculeType.RNA || moleculeType === MoleculeType.PNA
}

export function isNucleic(moleculeType: MoleculeType) {
    return moleculeType === MoleculeType.DNA || moleculeType === MoleculeType.RNA || moleculeType === MoleculeType.PNA
}

export function isProtein(moleculeType: MoleculeType) {
    return moleculeType === MoleculeType.protein
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
])

export type SecondaryStructureType = BitFlags<SecondaryStructureType.Flag>
export namespace SecondaryStructureType {
    export const is: (ss: SecondaryStructureType, f: Flag) => boolean = BitFlags.has
    export const create: (fs: Flag) => SecondaryStructureType = BitFlags.create

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
    }

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
    }

    export const SecondaryStructureStride: { [value: string]: number } = {
        H: Flag.Helix | Flag.HelixAlpha,  // Alpha helix
        G: Flag.Helix | Flag.Helix3Ten,  // 3-10 helix
        I: Flag.Helix | Flag.HelixPi,  // PI-helix
        E: Flag.Beta | Flag.BetaSheet,  // Extended conformation
        B: Flag.Beta | Flag.BetaStrand,  // Isolated bridge
        T: Flag.Turn,  // Turn
        C: Flag.NA,  // Coil (none of the above)
    }

    export const SecondaryStructureDssp: { [value: string]: number } = {
        H: Flag.Helix | Flag.HelixAlpha,  // alpha-helix
        B: Flag.Beta | Flag.BetaStrand,  // residue in isolated beta-bridge
        E: Flag.Beta | Flag.BetaSheet,  // extended strand, participates in beta ladder
        G: Flag.Helix | Flag.Helix3Ten,  // 3-helix (310 helix)
        I: Flag.Helix | Flag.HelixPi,  // 5 helix (pi-helix)
        T: Flag.Turn,  // hydrogen bonded turn
        S: Flag.Bend,  // bend
    }
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
}
export const DefaultMaxAsa = 121.0

export const VdwRadii = {
    'H': 1.1,
    'HE': 1.4,
    'LI': 1.81,
    'BE': 1.53,
    'B': 1.92,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'NE': 1.54,
    'NA': 2.27,
    'MG': 1.73,
    'AL': 1.84,
    'SI': 2.1,
    'P': 1.8,
    'S': 1.8,
    'CL': 1.75,
    'AR': 1.88,
    'K': 2.75,
    'CA': 2.31,
    'SC': 2.3,
    'TI': 2.15,
    'V': 2.05,
    'CR': 2.05,
    'MN': 2.05,
    'FE': 2.05,
    'CO': 2.0,
    'NI': 2.0,
    'CU': 2.0,
    'ZN': 2.1,
    'GA': 1.87,
    'GE': 2.11,
    'AS': 1.85,
    'SE': 1.9,
    'BR': 1.83,
    'KR': 2.02,
    'RB': 3.03,
    'SR': 2.49,
    'Y': 2.4,
    'ZR': 2.3,
    'NB': 2.15,
    'MO': 2.1,
    'TC': 2.05,
    'RU': 2.05,
    'RH': 2.0,
    'PD': 2.05,
    'AG': 2.1,
    'CD': 2.2,
    'IN': 2.2,
    'SN': 1.93,
    'SB': 2.17,
    'TE': 2.06,
    'I': 1.98,
    'XE': 2.16,
    'CS': 3.43,
    'BA': 2.68,
    'LA': 2.5,
    'CE': 2.48,
    'PR': 2.47,
    'ND': 2.45,
    'PM': 2.43,
    'SM': 2.42,
    'EU': 2.4,
    'GD': 2.38,
    'TB': 2.37,
    'DY': 2.35,
    'HO': 2.33,
    'ER': 2.32,
    'TM': 2.3,
    'YB': 2.28,
    'LU': 2.27,
    'HF': 2.25,
    'TA': 2.2,
    'W': 2.1,
    'RE': 2.05,
    'OS': 2.0,
    'IR': 2.0,
    'PT': 2.05,
    'AU': 2.1,
    'HG': 2.05,
    'TL': 1.96,
    'PB': 2.02,
    'BI': 2.07,
    'PO': 1.97,
    'AT': 2.02,
    'RN': 2.2,
    'FR': 3.48,
    'RA': 2.83,
    'AC': 2.0,
    'TH': 2.4,
    'PA': 2.0,
    'U': 2.3,
    'NP': 2.0,
    'PU': 2.0,
    'AM': 2.0,
    'CM': 2.0,
    'BK': 2.0,
    'CF': 2.0,
    'ES': 2.0,
    'FM': 2.0,
    'MD': 2.0,
    'NO': 2.0,
    'LR': 2.0,
    'RF': 2.0,
    'DB': 2.0,
    'SG': 2.0,
    'BH': 2.0,
    'HS': 2.0,
    'MT': 2.0,
    'DS': 2.0,
    'RG': 2.0,
    'CN': 2.0,
    'UUT': 2.0,
    'FL': 2.0,
    'UUP': 2.0,
    'LV': 2.0,
    'UUH': 2.0
}
export const DefaultVdwRadius = 2.0

export type LinkType = BitFlags<LinkType.Flag>
export namespace LinkType {
    export const is: (b: LinkType, f: Flag) => boolean = BitFlags.has
    export const enum Flag {
        None                 = 0x0,
        Covalent             = 0x1,
        MetallicCoordination = 0x2,
        Hydrogen             = 0x4,
        Ionic                = 0x8,
        Sulfide              = 0x10,
        Aromatic             = 0x20,
        Computed             = 0x40
        // currently at most 16 flags are supported!!
    }

    export function create(flags: Flag): LinkType {
        return BitFlags.create(flags);
    }

    export function isCovalent(flags: LinkType.Flag) {
        return (flags & LinkType.Flag.Covalent) !== 0;
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
    'ALA': [ 0.17, 0.50, 0.33 ],
    'ARG': [ 0.81, 1.81, 1.00 ],
    'ASN': [ 0.42, 0.85, 0.43 ],
    'ASP': [ 1.23, 3.64, 2.41 ],
    'ASH': [ -0.07, 0.43, 0.50 ],
    'CYS': [ -0.24, -0.02, 0.22 ],
    'GLN': [ 0.58, 0.77, 0.19 ],
    'GLU': [ 2.02, 3.63, 1.61 ],
    'GLH': [ -0.01, 0.11, 0.12 ],
    'GLY': [ 0.01, 1.15, 1.14 ],
    // "His+": [  0.96,  2.33,  1.37 ],
    'HIS': [ 0.17, 0.11, -0.06 ],
    'ILE': [ -0.31, -1.12, -0.81 ],
    'LEU': [ -0.56, -1.25, -0.69 ],
    'LYS': [ 0.99, 2.80, 1.81 ],
    'MET': [ -0.23, -0.67, -0.44 ],
    'PHE': [ -1.13, -1.71, -0.58 ],
    'PRO': [ 0.45, 0.14, -0.31 ],
    'SER': [ 0.13, 0.46, 0.33 ],
    'THR': [ 0.14, 0.25, 0.11 ],
    'TRP': [ -1.85, -2.09, -0.24 ],
    'TYR': [ -0.94, -0.71, 0.23 ],
    'VAL': [ 0.07, -0.46, -0.53 ]
  }
  export const DefaultResidueHydrophobicity = [ 0.00, 0.00, 0.00 ]