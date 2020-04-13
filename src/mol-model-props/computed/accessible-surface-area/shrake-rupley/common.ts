/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Structure } from '../../../../mol-model/structure';
import { Vec3 } from '../../../../mol-math/linear-algebra';

export interface ShrakeRupleyContext {
    structure: Structure,
    spherePoints: Vec3[],
    probeSize: number,
    nonPolymer: boolean,
    scalingConstant: number,
    maxLookupRadius: number,
    atomRadiusType: Int8Array,
    serialResidueIndex: Int32Array,
    /** Accessible surface area values */
    area: Float32Array
}

/** Chothia's amino acid and nucleotide atom vdw radii */
export const VdWLookup = [
    -1.0, // 0: missing
    1.76, // 1: trigonal C
    1.87, // 2: tetrahedral C
    1.65, // 3: trigonal N
    1.50, // 4: tetrahedral N
    1.40, // 5: O
    1.85, // 6: S
    1.80, // 7: C (nucleic)
    1.60, // 8: N (nucleic)
    1.40  // 9: P (nucleic)
]; // can still be appended on-the-fly for rare elements like selenium

/** Maximum accessible surface area observed for amino acids. Taken from: http://dx.doi.org/10.1371/journal.pone.0080635 */
export const MaxAsa: { [k: string]: number } = {
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