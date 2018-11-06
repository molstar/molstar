/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSphereVisual, ElementSphereParams } from '../visual/element-sphere';
import { UnitsRepresentation } from '../units-representation';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { StructureRepresentation } from '../index';
import { Representation } from 'mol-repr';

export const SpacefillParams = {
    ...ElementSphereParams
}
export const DefaultSpacefillProps = PD.paramDefaultValues(SpacefillParams)
export type SpacefillProps = typeof DefaultSpacefillProps

export type SpacefillRepresentation = StructureRepresentation<SpacefillProps>

export function SpacefillRepresentation(): SpacefillRepresentation {
    return Representation.createMulti('Spacefill', SpacefillParams, DefaultSpacefillProps, [
        UnitsRepresentation('Sphere mesh', ElementSphereVisual)
    ] as StructureRepresentation<SpacefillProps>[])
}