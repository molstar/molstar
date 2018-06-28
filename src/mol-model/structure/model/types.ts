/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import BitFlags from 'mol-util/bit-flags'

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
    'unknown', 'polymer', 'non-polymer', 'macrolide', 'water'
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
    /** sacharide, e.g. component type included in `SaccharideComponentTypeNames` */
    saccharide
}

/** Chemical component types as defined in the mmCIF CCD */
export enum ComponentType {
    // protein
    'D-peptide linking', 'L-peptide linking', 'D-peptide NH3 amino terminus',
    'L-peptide NH3 amino terminus', 'D-peptide COOH carboxy terminus',
    'L-peptide COOH carboxy terminus', 'peptide linking', 'peptide-like',
    'L-gamma-peptide, C-delta linking', 'D-gamma-peptide, C-delta linking',
    'L-beta-peptide, C-gamma linking', 'D-beta-peptide, C-gamma linking',

    // DNA
    'DNA linking', 'L-DNA linking', 'DNA OH 5 prime terminus', 'DNA OH 3 prime terminus',

    // RNA
    'RNA linking', 'L-RNA linking', 'RNA OH 5 prime terminus', 'RNA OH 3 prime terminus',

    // sacharide
    'D-saccharide 1,4 and 1,4 linking', 'L-saccharide 1,4 and 1,4 linking',
    'D-saccharide 1,4 and 1,6 linking', 'L-saccharide 1,4 and 1,6 linking', 'L-saccharide',
    'D-saccharide', 'saccharide',

    'non-polymer', 'other'
}

/** Chemical component type names for protein */
export const ProteinComponentTypeNames = [
    'D-peptide linking', 'L-peptide linking', 'D-peptide NH3 amino terminus',
    'L-peptide NH3 amino terminus', 'D-peptide COOH carboxy terminus',
    'L-peptide COOH carboxy terminus', 'peptide linking', 'peptide-like',
    'L-gamma-peptide, C-delta linking', 'D-gamma-peptide, C-delta linking',
    'L-beta-peptide, C-gamma linking', 'D-beta-peptide, C-gamma linking',
]

/** Chemical component type names for DNA */
export const DNAComponentTypeNames = [
    'DNA linking', 'L-DNA linking', 'DNA OH 5 prime terminus', 'DNA OH 3 prime terminus',
]

/** Chemical component type names for RNA */
export const RNAComponentTypeNames = [
    'RNA linking', 'L-RNA linking', 'RNA OH 5 prime terminus', 'RNA OH 3 prime terminus',
]

/** Chemical component type names for saccharide */
export const SaccharideComponentTypeNames = [
    'D-saccharide 1,4 and 1,4 linking', 'L-saccharide 1,4 and 1,4 linking',
    'D-saccharide 1,4 and 1,6 linking', 'L-saccharide 1,4 and 1,6 linking', 'L-saccharide',
    'D-saccharide', 'saccharide',
]

/** Common names for water molecules */
export const WaterNames = [
    'SOL', 'WAT', 'HOH', 'H2O', 'W', 'DOD', 'D3O', 'TIP3', 'TIP4', 'SPC'
]

/** get the molecule type from component type and id */
export function getMoleculeType(compType: string, compId: string) {
    if (ProteinComponentTypeNames.includes(compType)) {
        return MoleculeType.protein
    } else if (RNAComponentTypeNames.includes(compType)) {
        return MoleculeType.RNA
    } else if (DNAComponentTypeNames.includes(compType)) {
        return MoleculeType.DNA
    } else if (SaccharideComponentTypeNames.includes(compType)) {
        return MoleculeType.saccharide
    } else if (WaterNames.includes(compId)) {
        return MoleculeType.water
    } else if (IonNames.includes(compId)) {
        return MoleculeType.ion
    } else {
        return MoleculeType.unknown
    }
}

/**
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
export const IonNames = [
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
]

export interface SecondaryStructureType extends BitFlags<SecondaryStructureType.Flag> { }
export namespace SecondaryStructureType {
    export const Helix = ['h', 'g', 'i']
    export const Sheet = ['e', 'b']
    export const Turn = ['s', 't', 'l', '']

    export const is: (ss: SecondaryStructureType, f: Flag) => boolean = BitFlags.has
    export const create: (fs: Flag) => SecondaryStructureType = BitFlags.create

    export const enum Flag {
        None = 0x0,

        // category
        DoubleHelix = 0x1,
        Helix = 0x2,
        Beta = 0x4,
        Turn = 0x8,

        // category variant
        LeftHanded = 0x10,  // helix
        RightHanded = 0x20,

        ClassicTurn = 0x40,  // turn
        InverseTurn = 0x80,

        // sub-category
        HelixOther = 0x100,  // protein
        Helix27 = 0x200,
        Helix3Ten = 0x400,
        HelixAlpha = 0x800,
        HelixGamma = 0x1000,
        HelixOmega = 0x2000,
        HelixPi = 0x4000,
        HelixPolyproline = 0x8000,

        DoubleHelixOther = 0x10000,  // nucleic
        DoubleHelixZ = 0x20000,
        DoubleHelixA = 0x40000,
        DoubleHelixB = 0x80000,

        BetaOther = 0x100000,  // protein
        BetaStrand = 0x200000,  // single strand
        BetaSheet = 0x400000,  // multiple hydrogen bonded strands
        BetaBarell = 0x800000,  // closed series of sheets

        TurnOther = 0x1000000,  // protein
        Turn1 = 0x2000000,
        Turn2 = 0x4000000,
        Turn3 = 0x8000000,

        NA = 0x10000000,  // not applicable/available
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
        b: Flag.Beta | Flag.BetaStrand,  // Isolated bridge
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
        S: Flag.Turn,  // bend
    }
}

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

export interface LinkType extends BitFlags<LinkType.Flag> { }
export namespace LinkType {
    export const is: (b: LinkType, f: Flag) => boolean = BitFlags.has
    export const enum Flag {
        None                 = 0x0,
        Covalent             = 0x1,
        MetallicCoordination = 0x2,
        Hydrogen             = 0x4,
        Ion                  = 0x8,
        Sulfide              = 0x10,
        Aromatic             = 0x20,
        Computed             = 0x40
        // currently at most 16 flags are supported!!
    }

    export function isCovalent(flags: LinkType.Flag) {
        return (flags & LinkType.Flag.Covalent) !== 0;
    }
}