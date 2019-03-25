/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition'
import { BitFlags } from 'mol-util';

export interface AccessibleSurfaceArea {
    readonly atomRadius: ArrayLike<number>,
    readonly accessibleSurfaceArea: ArrayLike<number>,
    readonly relativeAccessibleSurfaceArea: ArrayLike<number>,
    readonly buried: Uint8Array
}

export const AccessibleSurfaceAreaComputationParams = {
    numberOfSpherePoints: PD.Numeric(92, {}, { description: 'number of sphere points to sample per atom: 92 (original paper), 960 (BioJava), 3000 (EPPIC) - see Shrake A, Rupley JA: Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol 1973.' }),
    probeSize: PD.Numeric(1.4, {}, { description: 'corresponds to the size of a water molecule: 1.4 (original paper), 1.5 (occassionally used)' }),
    buriedRasaThreshold: PD.Numeric(0.16, { min: 0.0, max: 1.0 }, { description: 'below this cutoff of relative accessible surface area a residue will be considered buried - see: Rost B, Sander C: Conservation and prediction of solvent accessibility in protein families. Proteins 1994.' }),
    nonPolymer: PD.Boolean(true, { description: 'Include non-polymer atoms in computation.' })
}

export namespace SolventAccessibility {
    export const is: (t: number, f: Flag) => boolean = BitFlags.has
    export const create: (f: Flag) => number = BitFlags.create
    export const enum Flag {
        _ = 0x0,
        BURIED = 0x1,
        ACCESSIBLE = 0x2
    }
}


export type AccessibleSurfaceAreaComputationParams = typeof AccessibleSurfaceAreaComputationParams