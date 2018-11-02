/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci } from 'mol-model/loci';
import { PolymerBackboneVisual, PolymerBackboneParams } from '../visual/polymer-backbone-cylinder';
import { getQualityProps } from '../../util';
import { paramDefaultValues } from 'mol-util/parameter';
import { UnitsRepresentation } from '../units-representation';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { StructureRepresentation } from '../index';

export const BackboneParams = {
    ...PolymerBackboneParams
}
export const DefaultBackboneProps = paramDefaultValues(BackboneParams)
export type BackboneProps = typeof DefaultBackboneProps

export type BackboneRepresentation = StructureRepresentation<BackboneProps>

export function BackboneRepresentation(): BackboneRepresentation {
    const traceRepr = UnitsRepresentation('Polymer backbone cylinder', PolymerBackboneVisual)

    let currentProps: BackboneProps
    return {
        label: 'Backbone',
        params: BackboneParams,
        get renderObjects() {
            return [ ...traceRepr.renderObjects ]
        },
        get props() {
            return { ...traceRepr.props }
        },
        createOrUpdate: (props: Partial<BackboneProps> = {}, structure?: Structure) => {
            const qualityProps = getQualityProps(Object.assign({}, currentProps, props), structure)
            currentProps = Object.assign({}, DefaultBackboneProps, currentProps, props, qualityProps)
            return Task.create('BackboneRepresentation', async ctx => {
                await traceRepr.createOrUpdate(currentProps, structure).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            return traceRepr.getLoci(pickingId)
        },
        mark: (loci: Loci, action: MarkerAction) => {
            return traceRepr.mark(loci, action)
        },
        destroy() {
            traceRepr.destroy()
        }
    }
}