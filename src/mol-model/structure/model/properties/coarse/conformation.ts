/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import UUID from 'mol-util/uuid'
import { Mat3 } from 'mol-math/linear-algebra';

export interface CoarseConformation {
    id: UUID,
    spheres: CoarseSphereConformation,
    gaussians: CoarseGaussianConformation
}

export interface CoarseSphereConformation {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,
    radius: ArrayLike<number>,
    rmsf: ArrayLike<number>
}

export interface CoarseGaussianConformation {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,
    weight: ArrayLike<number>,
    covariance_matrix: ArrayLike<Mat3>
}