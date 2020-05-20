/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Structure } from '../../../../mol-model/structure';
import { Vec3 } from '../../../../mol-math/linear-algebra';

export interface ANVILContext {
    structure: Structure,
    numberOfSpherePoints: number,
    stepSize: number,
    minThickness: number,
    maxThickness: number,
    afilter: number,
    membranePointDensity: number,

    offsets: ArrayLike<number>,
    exposed: ArrayLike<boolean>,
    centroid: Vec3,
    extent: number
};