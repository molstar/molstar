/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Structure } from '../../../../mol-model/structure';
import { Topology } from '../ANVIL';
import { Vec3 } from '../../../../mol-math/linear-algebra';

export interface ANVILContext {
    structure: Structure,
    numberOfSpherePoints: number,
    stepSize: number,
    minThickness: number,
    maxThickness: number,
    afilter: number,
    membranePointDensity: number,
    centerOfMass: Vec3,
    maxExtent: number,
    serialResidueIndex: Int32Array,
    exposure: ArrayLike<Topology>
}