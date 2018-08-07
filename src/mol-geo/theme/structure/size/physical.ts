/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Unit, StructureProperties } from 'mol-model/structure';
import { StructureSizeDataProps } from '.';
import { createAttributeSize } from '../../../util/size-data';

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
export function physicalSizeData(factor: number, props: StructureSizeDataProps) {
    const { group, vertexMap } = props
    const unit = group.units[0]
    const elements = group.elements;
    const radius = getPhysicalRadius(unit)
    const l = StructureElement.create()
    l.unit = unit
    return createAttributeSize({
        sizeFn: (elementIdx: number) => {
            l.element = elements[elementIdx]
            return radius(l) * factor
        },
        vertexMap
    })
}