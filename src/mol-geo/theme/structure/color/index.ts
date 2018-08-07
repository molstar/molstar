/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit } from 'mol-model/structure';

export interface StructureColorDataProps {
    group: Unit.SymmetryGroup,
    elementCount: number
}

export { elementIndexColorData } from './element-index'
export { chainIdElementColorData } from './chain-id'
export { elementSymbolColorData } from './element-symbol'
export { instanceIndexColorData } from './instance-index'