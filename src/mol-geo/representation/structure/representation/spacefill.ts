/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UnitsRepresentation } from '..';
import { ElementSphereVisual, DefaultElementSphereProps } from '../visual/element-sphere';
import { StructureRepresentation } from '../units-representation';

export const DefaultSpacefillProps = {
    ...DefaultElementSphereProps
}
export type SpacefillProps = typeof DefaultSpacefillProps

export type SpacefillRepresentation = StructureRepresentation<SpacefillProps>

export function SpacefillRepresentation(): SpacefillRepresentation {
    return UnitsRepresentation(ElementSphereVisual)
}