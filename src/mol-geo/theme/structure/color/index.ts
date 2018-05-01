/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementGroup, Unit } from 'mol-model/structure';
import VertexMap from '../../../shape/vertex-map';

export interface StructureColorDataProps {
    units: ReadonlyArray<Unit>,
    elementGroup: ElementGroup,
    vertexMap: VertexMap
}

export { elementIndexColorData } from './element-index'
export { chainIdColorData } from './chain-id'
export { elementSymbolColorData } from './element-symbol'
export { instanceIndexColorData } from './instance-index'