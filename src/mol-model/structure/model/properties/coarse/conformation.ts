/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import UUID from '../../../../../mol-util/uuid';
import { Mat3 } from '../../../../../mol-math/linear-algebra';

export interface CoarseConformation {
    id: UUID,
    spheres: CoarseSphereConformation,
    gaussians: CoarseGaussianConformation
}

export interface CoarseSphereConformation {
    /**
     * The x coordinate in angstroms specified according to a set of orthogonal Cartesian axes.
     */
    x: ArrayLike<number>,
    /**
     * The y coordinate in angstroms specified according to a set of orthogonal Cartesian axes.
     */
    y: ArrayLike<number>,
    /**
     * The z coordinate in angstroms specified according to a set of orthogonal Cartesian axes.
     */
    z: ArrayLike<number>,
    /**
     * The radius associated with the primitive sphere object at this position.
     */
    radius: ArrayLike<number>,
    /**
     * The Root Mean Square Fluctuation (RMSF) observed in the primitive sphere object at this position.
     */
    rmsf: ArrayLike<number>
}

export interface CoarseGaussianConformation {
    /**
     * The x coordinate in angstroms specified according to a set of orthogonal Cartesian axes.
     */
    x: ArrayLike<number>,
    /**
     * The y coordinate in angstroms specified according to a set of orthogonal Cartesian axes.
     */
    y: ArrayLike<number>,
    /**
     * The z coordinate in angstroms specified according to a set of orthogonal Cartesian axes.
     */
    z: ArrayLike<number>,
    /**
     * The weight of the gaussian object.
     */
    weight: ArrayLike<number>,
    /**
     * Data item of the covariance matrix representing the Gaussian object.
     */
    covariance_matrix: ArrayLike<Mat3>
}