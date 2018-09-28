/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UnitsVisual, VisualUpdateState } from '..';
import { createElementSphereMesh, markElement, getElementLoci, StructureElementIterator } from './util/element';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';

export const ElementSphereParams = {
    UnitsMe
}
export const DefaultElementSphereProps = {
    ...DefaultUnitsMeshProps,
    detail: 0
}
export type ElementSphereProps = typeof DefaultElementSphereProps

export function ElementSphereVisual(): UnitsVisual<ElementSphereProps> {
    return UnitsMeshVisual<ElementSphereProps>({
        defaultProps: DefaultElementSphereProps,
        createMesh: createElementSphereMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: ElementSphereProps, currentProps: ElementSphereProps) => {
            state.createGeometry = newProps.detail !== currentProps.detail
        }
    })
}