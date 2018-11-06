/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UnitsVisual } from '../index';
import { VisualUpdateState } from '../../util';
import { createElementSphereMesh, markElement, getElementLoci, StructureElementIterator } from './util/element';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';

export const ElementSphereParams = {
    ...UnitsMeshParams,
    sizeTheme: PD.Select<SizeThemeName>('Size Theme', '', 'physical', SizeThemeOptions),
    sizeValue: PD.Numeric('Size Value', '', 0.2, 0, 10, 0.1),
    sizeFactor: PD.Numeric('Size Factor', '', 1, 0, 10, 0.1),
    detail: PD.Numeric('Sphere Detail', '', 0, 0, 3, 1),
}
export const DefaultElementSphereProps = PD.getDefaultValues(ElementSphereParams)
export type ElementSphereProps = typeof DefaultElementSphereProps

export function ElementSphereVisual(): UnitsVisual<ElementSphereProps> {
    return UnitsMeshVisual<ElementSphereProps>({
        defaultProps: DefaultElementSphereProps,
        createGeometry: createElementSphereMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: ElementSphereProps, currentProps: ElementSphereProps) => {
            state.createGeometry = newProps.detail !== currentProps.detail
        }
    })
}