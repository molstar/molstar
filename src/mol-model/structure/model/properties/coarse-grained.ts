/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database as mmCIF } from 'mol-io/reader/cif/schema/mmcif'
import { Tensor } from 'mol-math/linear-algebra';
import { Column } from 'mol-data/db';

interface CoarseGrained {
    isDefined: boolean,
    modelList: mmCIF['ihm_model_list'],
    spheres: CoarseGrained.Spheres,
    gaussians: CoarseGrained.Gaussians
}

namespace CoarseGrained {
    export const Empty: CoarseGrained = { isDefined: false } as any;

    interface Site {
        // index to the Model.hierarchy.entities table
        entityKey: number,
        // index to the CoarseGrained.modelList table
        modelKey: number,

        asym_id: string,
        seq_id_begin: number,
        seq_id_end: number,
        x: number,
        y: number,
        z: number
    }

    export interface Sphere extends Site {
        radius: number,
        rmsf: number
    }

    export interface Gaussian extends Site {
        weight: number,
        covarianceMatrix: Tensor.Data
    }

    export type Spheres = { count: number} & { [P in keyof Sphere]: Column<Sphere[P]> }
    export type Gaussians = { count: number} & { [P in keyof Gaussian]: Column<Gaussian[P]> }
}

export default CoarseGrained;