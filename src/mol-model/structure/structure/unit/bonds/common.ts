/**
 * Copyright (c) 2017-2020 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSymbol } from '../../../model/types';

export interface BondComputationProps {
    forceCompute: boolean
    noCompute: boolean
}
export const DefaultBondComputationProps: BondComputationProps = {
    forceCompute: false,
    noCompute: false
};

// H,D,T are all mapped to H
const __ElementIndex: { [e: string]: number | undefined } = { 'H': 0, 'h': 0, 'D': 0, 'd': 0, 'T': 0, 't': 0, 'He': 2, 'HE': 2, 'he': 2, 'Li': 3, 'LI': 3, 'li': 3, 'Be': 4, 'BE': 4, 'be': 4, 'B': 5, 'b': 5, 'C': 6, 'c': 6, 'N': 7, 'n': 7, 'O': 8, 'o': 8, 'F': 9, 'f': 9, 'Ne': 10, 'NE': 10, 'ne': 10, 'Na': 11, 'NA': 11, 'na': 11, 'Mg': 12, 'MG': 12, 'mg': 12, 'Al': 13, 'AL': 13, 'al': 13, 'Si': 14, 'SI': 14, 'si': 14, 'P': 15, 'p': 15, 'S': 16, 's': 16, 'Cl': 17, 'CL': 17, 'cl': 17, 'Ar': 18, 'AR': 18, 'ar': 18, 'K': 19, 'k': 19, 'Ca': 20, 'CA': 20, 'ca': 20, 'Sc': 21, 'SC': 21, 'sc': 21, 'Ti': 22, 'TI': 22, 'ti': 22, 'V': 23, 'v': 23, 'Cr': 24, 'CR': 24, 'cr': 24, 'Mn': 25, 'MN': 25, 'mn': 25, 'Fe': 26, 'FE': 26, 'fe': 26, 'Co': 27, 'CO': 27, 'co': 27, 'Ni': 28, 'NI': 28, 'ni': 28, 'Cu': 29, 'CU': 29, 'cu': 29, 'Zn': 30, 'ZN': 30, 'zn': 30, 'Ga': 31, 'GA': 31, 'ga': 31, 'Ge': 32, 'GE': 32, 'ge': 32, 'As': 33, 'AS': 33, 'as': 33, 'Se': 34, 'SE': 34, 'se': 34, 'Br': 35, 'BR': 35, 'br': 35, 'Kr': 36, 'KR': 36, 'kr': 36, 'Rb': 37, 'RB': 37, 'rb': 37, 'Sr': 38, 'SR': 38, 'sr': 38, 'Y': 39, 'y': 39, 'Zr': 40, 'ZR': 40, 'zr': 40, 'Nb': 41, 'NB': 41, 'nb': 41, 'Mo': 42, 'MO': 42, 'mo': 42, 'Tc': 43, 'TC': 43, 'tc': 43, 'Ru': 44, 'RU': 44, 'ru': 44, 'Rh': 45, 'RH': 45, 'rh': 45, 'Pd': 46, 'PD': 46, 'pd': 46, 'Ag': 47, 'AG': 47, 'ag': 47, 'Cd': 48, 'CD': 48, 'cd': 48, 'In': 49, 'IN': 49, 'in': 49, 'Sn': 50, 'SN': 50, 'sn': 50, 'Sb': 51, 'SB': 51, 'sb': 51, 'Te': 52, 'TE': 52, 'te': 52, 'I': 53, 'i': 53, 'Xe': 54, 'XE': 54, 'xe': 54, 'Cs': 55, 'CS': 55, 'cs': 55, 'Ba': 56, 'BA': 56, 'ba': 56, 'La': 57, 'LA': 57, 'la': 57, 'Ce': 58, 'CE': 58, 'ce': 58, 'Pr': 59, 'PR': 59, 'pr': 59, 'Nd': 60, 'ND': 60, 'nd': 60, 'Pm': 61, 'PM': 61, 'pm': 61, 'Sm': 62, 'SM': 62, 'sm': 62, 'Eu': 63, 'EU': 63, 'eu': 63, 'Gd': 64, 'GD': 64, 'gd': 64, 'Tb': 65, 'TB': 65, 'tb': 65, 'Dy': 66, 'DY': 66, 'dy': 66, 'Ho': 67, 'HO': 67, 'ho': 67, 'Er': 68, 'ER': 68, 'er': 68, 'Tm': 69, 'TM': 69, 'tm': 69, 'Yb': 70, 'YB': 70, 'yb': 70, 'Lu': 71, 'LU': 71, 'lu': 71, 'Hf': 72, 'HF': 72, 'hf': 72, 'Ta': 73, 'TA': 73, 'ta': 73, 'W': 74, 'w': 74, 'Re': 75, 'RE': 75, 're': 75, 'Os': 76, 'OS': 76, 'os': 76, 'Ir': 77, 'IR': 77, 'ir': 77, 'Pt': 78, 'PT': 78, 'pt': 78, 'Au': 79, 'AU': 79, 'au': 79, 'Hg': 80, 'HG': 80, 'hg': 80, 'Tl': 81, 'TL': 81, 'tl': 81, 'Pb': 82, 'PB': 82, 'pb': 82, 'Bi': 83, 'BI': 83, 'bi': 83, 'Po': 84, 'PO': 84, 'po': 84, 'At': 85, 'AT': 85, 'at': 85, 'Rn': 86, 'RN': 86, 'rn': 86, 'Fr': 87, 'FR': 87, 'fr': 87, 'Ra': 88, 'RA': 88, 'ra': 88, 'Ac': 89, 'AC': 89, 'ac': 89, 'Th': 90, 'TH': 90, 'th': 90, 'Pa': 91, 'PA': 91, 'pa': 91, 'U': 92, 'u': 92, 'Np': 93, 'NP': 93, 'np': 93, 'Pu': 94, 'PU': 94, 'pu': 94, 'Am': 95, 'AM': 95, 'am': 95, 'Cm': 96, 'CM': 96, 'cm': 96, 'Bk': 97, 'BK': 97, 'bk': 97, 'Cf': 98, 'CF': 98, 'cf': 98, 'Es': 99, 'ES': 99, 'es': 99, 'Fm': 100, 'FM': 100, 'fm': 100, 'Md': 101, 'MD': 101, 'md': 101, 'No': 102, 'NO': 102, 'no': 102, 'Lr': 103, 'LR': 103, 'lr': 103, 'Rf': 104, 'RF': 104, 'rf': 104, 'Db': 105, 'DB': 105, 'db': 105, 'Sg': 106, 'SG': 106, 'sg': 106, 'Bh': 107, 'BH': 107, 'bh': 107, 'Hs': 108, 'HS': 108, 'hs': 108, 'Mt': 109, 'MT': 109, 'mt': 109 };

// Increased P (15) threshold from 1.9 to 2.0 (e.g. for G16 in 1O08)
const __ElementBondThresholds: { [e: number]: number | undefined } = { 0: 1.42, 1: 1.42, 3: 2.7, 4: 2.7, 6: 1.75, 7: 1.6, 8: 1.52, 11: 2.7, 12: 2.7, 13: 2.7, 14: 1.9, 15: 2.0, 16: 1.9, 17: 1.8, 19: 2.7, 20: 2.7, 21: 2.7, 22: 2.7, 23: 2.7, 24: 2.7, 25: 2.7, 26: 2.7, 27: 2.7, 28: 2.7, 29: 2.7, 30: 2.7, 31: 2.7, 33: 2.68, 37: 2.7, 38: 2.7, 39: 2.7, 40: 2.7, 41: 2.7, 42: 2.7, 43: 2.7, 44: 2.7, 45: 2.7, 46: 2.7, 47: 2.7, 48: 2.7, 49: 2.7, 50: 2.7, 55: 2.7, 56: 2.7, 57: 2.7, 58: 2.7, 59: 2.7, 60: 2.7, 61: 2.7, 62: 2.7, 63: 2.7, 64: 2.7, 65: 2.7, 66: 2.7, 67: 2.7, 68: 2.7, 69: 2.7, 70: 2.7, 71: 2.7, 72: 2.7, 73: 2.7, 74: 2.7, 75: 2.7, 76: 2.7, 77: 2.7, 78: 2.7, 79: 2.7, 80: 2.7, 81: 2.7, 82: 2.7, 83: 2.7, 87: 2.7, 88: 2.7, 89: 2.7, 90: 2.7, 91: 2.7, 92: 2.7, 93: 2.7, 94: 2.7, 95: 2.7, 96: 2.7, 97: 2.7, 98: 2.7, 99: 2.7, 100: 2.7, 101: 2.7, 102: 2.7, 103: 2.7, 104: 2.7, 105: 2.7, 106: 2.7, 107: 2.7, 108: 2.7, 109: 2.88 };

/**
 * Increased N-N (112) threshold from 1.55 to 1.6 (e.g. for 0QH in 1BMA)
 *
 * More experimentally observed bond length here (https://cccbdb.nist.gov/expbondlengths1x.asp)
 *
 * H-H (https://cccbdb.nist.gov/expbondlengths2x.asp?descript=rHH)
 * - 0.741
 *
 * Changed C-H (27) to 1.2
 * C-H (https://cccbdb.nist.gov/expbondlengths2x.asp?descript=rCH)
 * - Average    1.091 (+/- 0.017)
 * - Min        0.931
 * - Max        1.140
 *
 * Changed N-H (35) to 1.15
 * N-H (https://cccbdb.nist.gov/expbondlengths2x.asp?descript=rNH)
 * - Average    1.009 (+/- 0.043)
 * - Min        0.836
 * - Max        1.090
 *
 * Changed O-H (44) to 1.1
 * O-H (https://cccbdb.nist.gov/expbondlengths2x.asp?descript=rOH)
 * - Average    0.967 (+/- 0.022)
 * - Min        0.912
 * - Max        1.033
 *
 * Added P-H (135) as 1.47
 * P-H (https://cccbdb.nist.gov/expbondlengths2x.asp?descript=rPH)
 * - Average    1.423 ((+/- 0.007)
 * - Min        0.912
 * - Max        1.033
 *
 * Added S-H (152) as 1.45
 * S-H (https://cccbdb.nist.gov/expbondlengths2x.asp?descript=rSH)
 * - Average    1.345 (+/- 0.020)
 * - Min        1.322
 * - Max        1.400
 *
 * Added Si-Si (420) as 2.37
 * (https://cccbdb.nist.gov/expbondlengths2x.asp?descript=rSiSi)
 */
const __ElementPairThresholds: { [e: number]: number | undefined } = { 0: 0.8, 20: 1.31, 27: 1.2, 35: 1.15, 44: 1.1, 54: 1, 60: 1.84, 72: 1.88, 84: 1.75, 85: 1.56, 86: 1.76, 98: 1.6, 99: 1.68, 100: 1.63, 112: 1.6, 113: 1.59, 114: 1.36, 129: 1.45, 135: 1.47, 144: 1.6, 152: 1.45, 170: 1.4, 180: 1.55, 202: 2.4, 222: 2.24, 224: 1.91, 225: 1.98, 243: 2.02, 269: 2, 293: 1.9, 420: 2.37, 480: 2.3, 512: 2.3, 544: 2.3, 612: 2.1, 629: 1.54, 665: 1, 813: 2.6, 854: 2.27, 894: 1.93, 896: 2.1, 937: 2.05, 938: 2.06, 981: 1.62, 1258: 2.68, 1309: 2.33, 1484: 1, 1763: 2.14, 1823: 2.48, 1882: 2.1, 1944: 1.72, 2380: 2.34, 3367: 2.44, 3733: 2.11, 3819: 2.6, 3821: 2.36, 4736: 2.75, 5724: 2.73, 5959: 2.63, 6519: 2.84, 6750: 2.87, 8991: 2.81 };

const __DefaultBondingRadius = 2.001;

export const MetalsSet = (function () {
    const metals = ['LI', 'NA', 'K', 'RB', 'CS', 'FR', 'BE', 'MG', 'CA', 'SR', 'BA', 'RA', 'AL', 'GA', 'IN', 'SN', 'TL', 'PB', 'BI', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'LA', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'AC', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR'];
    const set = new Set<number>();
    for (const m of metals) {
        set.add(__ElementIndex[m]!);
    }
    return set;
})();

function pair(a: number, b: number) {
    if (a < b) return (a + b) * (a + b + 1) / 2 + b;
    else return (a + b) * (a + b + 1) / 2 + a;
}

export function getElementIdx(e: ElementSymbol) {
    const i = __ElementIndex[e as any as string];
    if (i === void 0) return -1;
    return i;
}

export function getElementPairThreshold(i: number, j: number) {
    if (i < 0 || j < 0) return -1;
    const r = __ElementPairThresholds[pair(i, j)];
    if (r === void 0) return -1;
    return r;
}

export function getElementThreshold(i: number) {
    if (i < 0) return __DefaultBondingRadius;
    const r = __ElementBondThresholds[i];
    if (r === void 0) return __DefaultBondingRadius;
    return r;
}

const H_ID = __ElementIndex['H']!;
export function isHydrogen(i: number) {
    return i === H_ID;
}