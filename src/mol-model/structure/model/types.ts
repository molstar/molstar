/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
export interface ElementSymbol extends String { '@type': 'element-symbol' }
export function ElementSymbol(s: string): ElementSymbol {
    return _esCache[s] || s.toUpperCase();
}

export const enum EntityType {
    Unknown = 'unknown',
    Polymer = 'polymer',
    NonPolymer = 'non-polymer',
    Macrolide = 'macrolide',
    Water = 'water'
}

export const enum MoleculeType {
    Unknown,
    Water,
    Ion,
    Protein,
    RNA,
    DNA,
    Saccharide
}

export const enum BackboneType {
    Unknown,
    Protein,
    RNA,
    DNA,
    CgProtein,
    CgRNA,
    CgDNA
}

const _chemCompNonPolymer = ['NON-POLYMER'];
const _chemCompOther = ['OTHER'];
const _chemCompSaccharide = [
    'D-SACCHARIDE', 'D-SACCHARIDE 1,4 AND 1,4 LINKING', 'D-SACCHARIDE 1,4 AND 1,6 LINKING',
    'L-SACCHARIDE', 'L-SACCHARIDE 1,4 AND 1,4 LINKING', 'L-SACCHARIDE 1,4 AND 1,6 LINKING',
    'SACCHARIDE'
];

export const ChemComp = {
    Protein: [
        'D-BETA-PEPTIDE, C-GAMMA LINKING', 'D-GAMMA-PEPTIDE, C-DELTA LINKING',
        'D-PEPTIDE COOH CARBOXY TERMINUS', 'D-PEPTIDE NH3 AMINO TERMINUS', 'D-PEPTIDE LINKING',
        'L-BETA-PEPTIDE, C-GAMMA LINKING', 'L-GAMMA-PEPTIDE, C-DELTA LINKING',
        'L-PEPTIDE COOH CARBOXY TERMINUS', 'L-PEPTIDE NH3 AMINO TERMINUS', 'L-PEPTIDE LINKING',
        'PEPTIDE LINKING', 'PEPTIDE-LIKE'
    ],
    RNA: [
        'RNA OH 3 PRIME TERMINUS', 'RNA OH 5 PRIME TERMINUS', 'RNA LINKING'
    ],
    DNA: [
        'DNA OH 3 PRIME TERMINUS', 'DNA OH 5 PRIME TERMINUS', 'DNA LINKING',
        'L-DNA LINKING', 'L-RNA LINKING'
    ],
    Saccharide: _chemCompSaccharide,
    Other: _chemCompOther,
    NonPolymer: _chemCompNonPolymer,
    Hetero: _chemCompNonPolymer.concat(_chemCompOther, _chemCompSaccharide)
}

export interface SecondaryStructureType extends BitFlags<SecondaryStructureType.Flag> { }
export namespace SecondaryStructureType {
    export const Helix = ['h', 'g', 'i']
    export const Sheet = ['e', 'b']
    export const Turn = ['s', 't', 'l', '']

    export const is: (ss: SecondaryStructureType, f: Flag) => boolean = BitFlags.has

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