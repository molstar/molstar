/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database as mmCIF } from 'mol-io/reader/cif/schema/mmcif'
import { Tensor } from 'mol-math/linear-algebra';
import { Column } from 'mol-data/db';

interface CoarseGrainedHierarchy {
    isDefined: boolean,
    modelList: mmCIF['ihm_model_list'],
    spheres: CoarseGrainedHierarchy.Spheres,
    gaussians: CoarseGrainedHierarchy.Gaussians
}

namespace CoarseGrainedHierarchy {
    export const Empty: CoarseGrainedHierarchy = { isDefined: false } as any;

    export const enum ElementType { Sphere, Gaussian }

    export interface SiteBase {
        asym_id: string,
        seq_id_begin: number,
        seq_id_end: number
    }

    export interface Sphere extends SiteBase {
        radius: number,
        rmsf: number
    }

    export interface Gaussian extends SiteBase {
        weight: number,
        covariance_matrix: Tensor.Data
    }

    type Common = {
        count: number,
        x: ArrayLike<number>,
        y: ArrayLike<number>,
        z: ArrayLike<number>,
        modelKey: ArrayLike<number>,
        entityKey: ArrayLike<number>
    }

    export type SitesBase =  Common & { [P in keyof SiteBase]: Column<SiteBase[P]> }
    export type Spheres =  Common & { [P in keyof Sphere]: Column<Sphere[P]> }
    export type Gaussians = Common & { matrix_space: Tensor.Space } & { [P in keyof Gaussian]: Column<Gaussian[P]> }
}

export { CoarseGrainedHierarchy }