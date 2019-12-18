/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from './structure'
import { Bond } from './structure/structure/unit/bonds'
import { ShapeGroup } from './shape/shape';
import { Interactions } from '../mol-model-props/computed/interactions/interactions';

/** A null value Location */
export const NullLocation = { kind: 'null-location' as 'null-location' }
export type NullLocation = typeof NullLocation
export function isNullLocation(x: any): x is NullLocation {
    return !!x && x.kind === 'null-location';
}

export type Location = StructureElement.Location | Bond.Location | Interactions.Location | ShapeGroup.Location | NullLocation