/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UnitsVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { createElementSphereMesh, markElement, getElementLoci, StructureElementIterator } from './util/element';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { BuiltInSizeThemeName, BuiltInSizeThemeOptions } from 'mol-theme/size';

export const ElementSphereParams = {
    ...UnitsMeshParams,
    sizeTheme: PD.Select<BuiltInSizeThemeName>('Size Theme', '', 'physical', BuiltInSizeThemeOptions),
    sizeFactor: PD.Numeric('Size Factor', '', 1, 0, 10, 0.1),
    detail: PD.Numeric('Sphere Detail', '', 0, 0, 3, 1),
}
export type ElementSphereParams = typeof ElementSphereParams

export function ElementSphereVisual(): UnitsVisual<ElementSphereParams> {
    return UnitsMeshVisual<ElementSphereParams>({
        defaultProps: PD.getDefaultValues(ElementSphereParams),
        createGeometry: createElementSphereMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<ElementSphereParams>, currentProps: PD.Values<ElementSphereParams>) => {
            state.createGeometry = newProps.detail !== currentProps.detail
        }
    })
}