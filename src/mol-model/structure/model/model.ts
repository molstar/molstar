/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import UUID from 'mol-util/uuid'
import GridLookup from 'mol-math/geometry/grid-lookup'
import Format from './format'
import Hierarchy from './properties/hierarchy'
import Conformation from './properties/conformation'
import Symmetry from './properties/symmetry'
import Bonds from './properties/bonds'
import CoarseGrained from './properties/coarse-grained'

import computeBonds from './utils/compute-bonds'

import from_gro from './formats/gro'
import from_mmCIF from './formats/mmcif'

import { Annotations } from '../../annotations'

/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
interface Model extends Readonly<{
    id: UUID,

    modelNum: number,

    sourceData: Format,

    hierarchy: Hierarchy,
    conformation: Conformation,
    symmetry: Symmetry,
    coarseGrained: CoarseGrained,
    annotations: Annotations,

    atomCount: number,
}> {
    '@spatialLookup'?: GridLookup,
    '@bonds'?: Bonds
} { }

namespace Model {
    export function create(format: Format) {
        switch (format.kind) {
            case 'gro': return from_gro(format);
            case 'mmCIF': return from_mmCIF(format);
        }
    }
    export function spatialLookup(model: Model): GridLookup {
        if (model['@spatialLookup']) return model['@spatialLookup']!;
        const lookup = GridLookup(model.conformation);
        model['@spatialLookup'] = lookup;
        return lookup;
    }
    export function bonds(model: Model): Bonds {
        if (model['@bonds']) return model['@bonds']!;
        const bonds = computeBonds(model);
        model['@bonds'] = bonds;
        return bonds;
    }
}

export default Model