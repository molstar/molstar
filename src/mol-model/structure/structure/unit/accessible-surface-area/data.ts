/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition'

export interface AccessibleSurfaceArea {
    readonly atomRadius: ArrayLike<number>,
    readonly accessibleSurfaceArea: ArrayLike<number>,
    readonly relativeAccessibleSurfaceArea: ArrayLike<number>,
    readonly buried: any
}

export const AccessibleSurfaceAreaComputationParams = {
    numberOfSpherePoints: PD.Numeric(92, {}, { description: 'number of sphere points to sample per atom: 92 (original paper), 960 (BioJava), 3000 (EPPIC)' }),
    probeSize: PD.Numeric(1.4, {}, { description: 'corresponds to the size of a water molecule: 1.4 (original paper), 1.5 (occassionally used)' })
}
export type AccessibleSurfaceAreaComputationParams = typeof AccessibleSurfaceAreaComputationParams