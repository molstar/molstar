/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from 'mol-model/structure';
import { AssemblySymmetry } from 'mol-model-props/rcsb/symmetry';

export function RCSB_assemblySymmetry(model: Model) {
    return AssemblySymmetry.attachFromCifOrAPI(model)
}