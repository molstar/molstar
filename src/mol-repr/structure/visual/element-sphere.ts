/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UnitsVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { createElementSphereMesh, markElement, getElementLoci, StructureElementIterator, createElementSphereImpostor } from './util/element';
import { UnitsMeshVisual, UnitsMeshParams, UnitsSpheresVisual, UnitsSpheresParams } from '../units-visual';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { WebGLContext } from 'mol-gl/webgl/context';

export const ElementSphereParams = {
    ...UnitsMeshParams,
    ...UnitsSpheresParams,
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }),
}
export type ElementSphereParams = typeof ElementSphereParams

export function getElementSphereVisual(webgl?: WebGLContext) {
    return webgl && webgl.extensions.fragDepth ? ElementSphereImpostorVisual : ElementSphereMeshVisual
}

export function ElementSphereImpostorVisual(): UnitsVisual<ElementSphereParams> {
    return UnitsSpheresVisual<ElementSphereParams>({
        defaultProps: PD.getDefaultValues(ElementSphereParams),
        createGeometry: createElementSphereImpostor,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<ElementSphereParams>, currentProps: PD.Values<ElementSphereParams>) => {

        }
    })
}

export function ElementSphereMeshVisual(): UnitsVisual<ElementSphereParams> {
    return UnitsMeshVisual<ElementSphereParams>({
        defaultProps: PD.getDefaultValues(ElementSphereParams),
        createGeometry: createElementSphereMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<ElementSphereParams>, currentProps: PD.Values<ElementSphereParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail
            )
        }
    })
}