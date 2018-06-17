/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureRepresentation } from '.';
import { ElementSphereVisual, DefaultElementSphereProps } from './visual/element-sphere';

export const DefaultSpacefillProps = {
    ...DefaultElementSphereProps,
}
export type SpacefillProps = Partial<typeof DefaultSpacefillProps>

export function SpacefillRepresentation() {
    return StructureRepresentation(ElementSphereVisual)
}