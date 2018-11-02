/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PolymerBackboneVisual, PolymerBackboneParams } from '../visual/polymer-backbone-cylinder';
import { paramDefaultValues } from 'mol-util/parameter';
import { UnitsRepresentation } from '../units-representation';
import { StructureRepresentation } from '../index';
import { Representation } from 'mol-repr';

export const BackboneParams = {
    ...PolymerBackboneParams
}
export const DefaultBackboneProps = paramDefaultValues(BackboneParams)
export type BackboneProps = typeof DefaultBackboneProps

export type BackboneRepresentation = StructureRepresentation<BackboneProps>

export function BackboneRepresentation(): BackboneRepresentation {
    return Representation.createMulti('Backbone', BackboneParams, DefaultBackboneProps, [
        UnitsRepresentation('Polymer backbone cylinder', PolymerBackboneVisual)
    ] as StructureRepresentation<BackboneProps>[])
}