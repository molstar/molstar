/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementPointVisual, ElementPointParams } from '../visual/element-point';
import { UnitsRepresentation } from '../units-representation';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { StructureRepresentation } from '../index';
import { Representation } from 'mol-repr';

export const PointParams = {
    ...ElementPointParams,
}
export const DefaultPointProps = PD.paramDefaultValues(PointParams)
export type PointProps = typeof DefaultPointProps

export type PointRepresentation = StructureRepresentation<PointProps>

export function PointRepresentation(): PointRepresentation {
    return Representation.createMulti('Point', PointParams, DefaultPointProps, [
        UnitsRepresentation('Point', ElementPointVisual)
    ] as StructureRepresentation<PointProps>[])
}