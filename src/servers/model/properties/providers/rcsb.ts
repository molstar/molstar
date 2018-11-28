/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssemblySymmetry } from 'mol-model-props/rcsb/assembly-symmetry';
import { AttachModelProperty } from '../../property-provider';

export const RCSB_assemblySymmetry: AttachModelProperty = ({ model }) => {
    return AssemblySymmetry.attachFromCifOrAPI(model)
}