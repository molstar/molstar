/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Queries, Element } from 'mol-model/structure';

import { StructureColorDataProps } from '.';
import { createAttributeOrElementColor, ColorData } from '../../../util/color-data';
import { ColorScale } from 'mol-util/color';
import { Column } from 'mol-data/db';

function createChainIdMap(unit: Unit) {
    const map = new Map<string, number>()
    let index = 0

    let count: number
    let asym_id: Column<string>
    if (Unit.isAtomic(unit)) {
        asym_id = unit.model.atomicHierarchy.chains.label_asym_id
        count = unit.model.atomicHierarchy.chains._rowCount
    } else if (Unit.isCoarse(unit)) {
        asym_id = unit.coarseElements.asym_id
        count = unit.coarseElements.count
    } else {
        console.warn('Unknown unit type')
        return { map, count: index }
    }

    for (let i = 0; i < count; ++i) {
        const chainId = asym_id.value(i)
        if (map.get(chainId) === undefined) {
            map.set(chainId, index)
            index += 1
        }
    }
    return { map, count: index }
}

export function chainIdColorData(props: StructureColorDataProps, colorData?: ColorData) {
    const { group: { units, elements }, vertexMap } = props
    const unit = units[0]

    const { map, count } = createChainIdMap(unit)

    const domain = [ 0, count - 1 ]
    const scale = ColorScale.create({ domain })

    let asym_id: Element.Property<string>
    if (Unit.isAtomic(unit)) {
        asym_id = Queries.props.chain.label_asym_id
    } else if (Unit.isCoarse(unit)) {
        asym_id = Queries.props.coarse.asym_id
    }

    const l = Element.Location()
    l.unit = unit

    return createAttributeOrElementColor(vertexMap, {
        colorFn: (elementIdx: number) => {
            l.element = elements[elementIdx]
            return scale.color(map.get(asym_id(l)) || 0)
        },
        vertexMap
    }, colorData)
}