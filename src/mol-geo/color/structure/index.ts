/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementGroup, Unit } from 'mol-model/structure';
import { Mesh } from '../../shape/mesh';

export interface StructureColorDataProps {
    units: ReadonlyArray<Unit>,
    elementGroup: ElementGroup,
    mesh: Mesh
}