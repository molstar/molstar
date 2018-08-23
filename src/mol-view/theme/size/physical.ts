/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Unit, Link, ElementIndex } from 'mol-model/structure';
import { Location } from 'mol-model/location';
import { SizeThemeProps, SizeTheme } from '../size';
import { VdwRadius } from 'mol-model/structure/model/properties/atomic';

const DefaultSize = 1
const DefaultFactor = 1

export function getPhysicalRadius(unit: Unit, element: ElementIndex): number {
    if (Unit.isAtomic(unit)) {
        return VdwRadius(unit.model.atomicHierarchy.atoms.type_symbol.value(element))
    } else if (Unit.isSpheres(unit)) {
        return unit.model.coarseConformation.spheres.radius[element]
    } else {
        return 0
    }
}

/**
 * Create attribute data with the physical size of an element,
 * i.e. vdw for atoms and radius for coarse spheres
 */
export function PhysicalSizeTheme(props: SizeThemeProps): SizeTheme {
    const factor = props.factor || DefaultFactor

    function sizeFn(location: Location): number {
        let size: number
        if (StructureElement.isLocation(location)) {
            size = getPhysicalRadius(location.unit, location.element)
        } else if (Link.isLocation(location)) {
            size = getPhysicalRadius(location.aUnit, location.aUnit.elements[location.aIndex])
        } else {
            size = DefaultSize
        }
        return factor * size
    }

    return {
        kind: 'group',
        size: sizeFn
    }
}