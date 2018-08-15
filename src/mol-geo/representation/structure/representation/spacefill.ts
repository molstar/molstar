/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UnitsRepresentation } from '..';
import { ElementSphereVisual, DefaultElementSphereProps } from '../visual/element-sphere';

export const DefaultSpacefillProps = {
    ...DefaultElementSphereProps
}
export type SpacefillProps = typeof DefaultSpacefillProps

export function SpacefillRepresentation() {
    return UnitsRepresentation(ElementSphereVisual)
}