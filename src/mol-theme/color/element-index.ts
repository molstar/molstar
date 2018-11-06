/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale, Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { StructureElement, Link } from 'mol-model/structure';
import { OrderedSet } from 'mol-data/int';
import { ColorThemeProps, ColorTheme, LocationColor } from '../color';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives every element (atom or coarse sphere/gaussian) a unique color based on the position (index) of the element in the list of elements in the structure.'

export function ElementIndexColorTheme(props: ColorThemeProps): ColorTheme {
    let color: LocationColor
    let scale = ColorScale.create({ list: props.list, minLabel: 'Start', maxLabel: 'End' })

    if (props.structure) {
        const { units } = props.structure
        const unitCount = units.length
        const cummulativeElementCount = new Map<number, number>()
        const unitIdIndex = new Map<number, number>()

        let elementCount = 0
        for (let i = 0; i < unitCount; ++i) {
            cummulativeElementCount.set(i, elementCount)
            elementCount += units[i].elements.length
            unitIdIndex.set(units[i].id, i)
        }
        scale.setDomain(0, elementCount - 1)
        const scaleColor = scale.color

        color = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                const unitIndex = unitIdIndex.get(location.unit.id)!
                const unitElementIndex = OrderedSet.findPredecessorIndex(location.unit.elements, location.element)
                return scaleColor(cummulativeElementCount.get(unitIndex)! + unitElementIndex)
            } else if (Link.isLocation(location)) {
                const unitIndex = unitIdIndex.get(location.aUnit.id)!
                return scaleColor(cummulativeElementCount.get(unitIndex)! + location.aIndex)
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        features: { structure: true, list: true },
        granularity: 'groupInstance',
        color,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}