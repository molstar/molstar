/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Unit, StructureProperties, Link } from 'mol-model/structure';
import { Location } from 'mol-model/location';
import { SizeThemeProps, SizeTheme } from '../size';

const DefaultSize = 1
const DefaultFactor = 1

export function getPhysicalRadius(unit: Unit): StructureElement.Property<number> {
    if (Unit.isAtomic(unit)) {
        return StructureProperties.atom.vdw_radius
    } else if (Unit.isSpheres(unit)) {
        return StructureProperties.coarse.sphere_radius
    } else {
        return () => 0
    }
}

/**
 * Create attribute data with the physical size of an element,
 * i.e. vdw for atoms and radius for coarse spheres
 */
export function PhysicalSizeTheme(props: SizeThemeProps): SizeTheme {
    const factor = props.factor || DefaultFactor
    const l = StructureElement.create()

    function sizeFn(location: Location): number {
        if (StructureElement.isLocation(location)) {
            if (Unit.isAtomic(location.unit)) {
                const radius = getPhysicalRadius(location.unit)
                return factor * radius(location)
            }
        } else if (Link.isLocation(location)) {
            if (Unit.isAtomic(location.aUnit)) {
                const radius = getPhysicalRadius(location.aUnit)
                l.unit = location.aUnit
                l.element = location.aUnit.elements[location.aIndex]
                return factor * radius(l)
            }
        }
        return DefaultSize
    }

    return {
        kind: 'element',
        size: sizeFn
    }
}