/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type AminoAlphabet =
    | 'H' | 'R' | 'K' | 'I' | 'F' | 'L' | 'W' | 'A' | 'M' | 'P' | 'C' | 'N' | 'V' | 'G' | 'S' | 'Q' | 'Y' | 'D' | 'E' | 'T' | 'U' | 'O'
    | 'X' /** = Unknown */
    | '-' /** = Gap */

export type NuclecicAlphabet =
    | 'A' | 'C' | 'G' | 'T' | 'U'
    | 'X' /** = Unknown */
    | '-' /** = Gap */

// from NGL
const ProteinOneLetterCodes: { [name: string]: AminoAlphabet }  = {
    'HIS': 'H',
    'ARG': 'R',
    'LYS': 'K',
    'ILE': 'I',
    'PHE': 'F',
    'LEU': 'L',
    'TRP': 'W',
    'ALA': 'A',
    'MET': 'M',
    'PRO': 'P',
    'CYS': 'C',
    'ASN': 'N',
    'VAL': 'V',
    'GLY': 'G',
    'SER': 'S',
    'GLN': 'Q',
    'TYR': 'Y',
    'ASP': 'D',
    'GLU': 'E',
    'THR': 'T',

    'SEC': 'U',  // as per IUPAC definition
    'PYL': 'O',  // as per IUPAC definition
};

const DnaOneLetterCodes: { [name: string]: NuclecicAlphabet } = {
    'DA': 'A',
    'DC': 'C',
    'DG': 'G',
    'DT': 'T',
    'DU': 'U'
};

const RnaOneLetterCodes: { [name: string]: NuclecicAlphabet } = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'U'
};

export function getProteinOneLetterCode(residueName: string): AminoAlphabet {
    const code = ProteinOneLetterCodes[residueName];
    return code || 'X';
}

export function getRnaOneLetterCode(residueName: string): NuclecicAlphabet {
    const code = RnaOneLetterCodes[residueName];
    return code || 'X';
}

export function getDnaOneLetterCode(residueName: string): NuclecicAlphabet {
    const code = DnaOneLetterCodes[residueName];
    return code || 'X';
}