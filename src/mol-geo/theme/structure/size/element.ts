/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Element, Unit, Queries } from 'mol-model/structure';
import { StructureSizeDataProps } from '.';
import { createAttributeSize } from '../../../util/size-data';

/** Create attribute data with the size of an element, i.e. vdw for atoms and radius for coarse spheres */
export function elementSizeData(props: StructureSizeDataProps) {
    const { group, vertexMap } = props
    const unit = group.units[0]
    const elements = group.elements;
    let radius: Element.Property<number>
    if (Unit.isAtomic(unit)) {
        radius = Queries.props.atom.vdw_radius
    } else if (Unit.isSpheres(unit)) {
        radius = Queries.props.coarse.sphere_radius
    }
    const l = Element.Location()
    l.unit = unit
    return createAttributeSize({
        sizeFn: (elementIdx: number) => {
            l.element = elements[elementIdx]
            return radius(l)
        },
        vertexMap
    })
}