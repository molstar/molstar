/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UnitsVisual, VisualUpdateState } from '..';
import { createElementSphereMesh, markElement, getElementLoci, StructureElementIterator } from './util/element';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { NumberParam, paramDefaultValues, SelectParam } from 'mol-view/parameter';
import { SizeThemeName, SizeThemeOptions } from 'mol-view/theme/size';

export const ElementSphereParams = {
    ...UnitsMeshParams,
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'physical', SizeThemeOptions),
    sizeValue: NumberParam('Size Value', '', 1, 0, 20, 0.1),
    detail: NumberParam('Sphere Detail', '', 0, 0, 3, 1),
}
export const DefaultElementSphereProps = paramDefaultValues(ElementSphereParams)
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