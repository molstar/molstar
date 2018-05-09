/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureSymmetry } from 'mol-model/structure';
import VertexMap from '../../../shape/vertex-map';

export interface StructureSizeDataProps {
    group: StructureSymmetry.UnitGroup,
    vertexMap: VertexMap
}

export { elementSizeData } from './element'