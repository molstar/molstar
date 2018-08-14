/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UnitsVisual, MeshUpdateState } from '..';
import { createElementSphereMesh, markElement, getElementLoci } from './util/element';
import { StructureElementIterator } from './util/location-iterator';
import { UnitsMeshVisual, DefaultUnitsMeshProps } from '../units-visual';

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
        setUpdateState: (state: MeshUpdateState, newProps: ElementSphereProps, currentProps: ElementSphereProps) => {
            state.createMesh = newProps.detail !== currentProps.detail
        }
    })
}